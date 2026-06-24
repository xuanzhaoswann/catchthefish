import wolf

class extract_sam(wolf.Task):
    name = "extract_sam"
    inputs = {
        "id":None,
        "gs_clean_bam":None,
        "gs_clean_bai":None,
        "OncoBed":None,
        "IgBed":None,
        "min_map_quality_onco":30,
        "min_map_quality_ig":0,
        "exclude_flag_split":1540,
        "exclude_flag_mate":1548,
        "include_flag_mate":1,
        "add_gs_suffix":True
    }
    overrides = {"gs_clean_bam":"string", "gs_clean_bai":"string"}
    script = """
    set -euxo pipefail

    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

    # in case we need to reprocess a bit the URI
    bam="${gs_clean_bam}"
    bai="${gs_clean_bai}"

    cat ${OncoBed} | \
        awk '{ print "$7==\\""$1"\\"","&& $8>"$2,"&& $8<="$3}' | \
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ || /g' | \
        cat <(echo -n "{ if (") - | \
        awk '{print $0" ) { print $0 ;}}"}' > \
        code_filter_sam_MP_Onco.awk

    cat ${OncoBed} | \
        awk '{ print "$12==\\""$1"\\"","&& $13>"$2,"&& $13<="$3}' | \
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ || /g' | \
        cat <(echo -n "{ if (") - | \
        awk '{print $0" ) { print $0 ;}}"}' > \
        code_filter_sam_SA_Onco.awk

    cat ${IgBed} | \
        awk '{ print "$7==\\""$1"\\"","&& $8>"$2,"&& $8<="$3}' | \
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ || /g' | \
        cat <(echo -n "{ if (") - | \
        awk '{print $0" ) { print $0 ;}}"}' > \
        code_filter_sam_MP_Ig.awk

    cat ${IgBed} | \
        awk '{ print "$12==\\""$1"\\"","&& $13>"$2,"&& $13<="$3}' | \
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ || /g' | \
        cat <(echo -n "{ if (") - | \
        awk '{print $0" ) { print $0 ;}}"}' > \
        code_filter_sam_SA_Ig.awk

    samtools view \
        -q ${min_map_quality_ig} \
        -F ${exclude_flag_mate} \
        -f ${include_flag_mate} \
        -M -L ${IgBed} -X ${bam} ${bai} | \
        cut -f1-11 | \
        awk -f code_filter_sam_MP_Onco.awk | \
        awk '{print "Paired-Read",$0}' OFS='\t' \
        > ${id}_MP_reads_Ig_Onco.txt

    samtools view \
        -q ${min_map_quality_ig} \
        -F ${exclude_flag_split} \
        -M -L ${IgBed} -e '[SA]' -X ${bam} ${bai} | \
        cut -f1-12 | \
        sed 's/SA:Z://g' | \
        awk -F'\t' 'gsub(/,/,"\t",$12)' OFS='\t' | \
        awk -f code_filter_sam_SA_Onco.awk | \
        cut -f1-11 | \
        awk '{print "Split-Read",$0}' OFS='\t' \
        > ${id}_SA_reads_Ig_Onco.txt

    samtools view \
        -q ${min_map_quality_onco} \
        -F ${exclude_flag_mate} \
        -f ${include_flag_mate} \
        -M -L ${OncoBed} -X ${bam} ${bai} | \
        cut -f1-11 | \
        awk -f code_filter_sam_MP_Ig.awk | \
        awk '{print "Paired-Read",$0}' OFS='\t' \
        > ${id}_MP_reads_Onco_Ig.txt

    samtools view \
        -q ${min_map_quality_onco} \
        -F ${exclude_flag_split} \
        -M -L ${OncoBed} -e '[SA]' -X ${bam} ${bai} | \
        cut -f1-12 | \
        sed 's/SA:Z://g' | \
        awk -F'\t' 'gsub(/,/,"\t",$12)' OFS='\t' | \
        awk -f code_filter_sam_SA_Ig.awk | \
        cut -f1-11 | \
        awk '{print "Split-Read",$0}' OFS='\t' \
        > ${id}_SA_reads_Onco_Ig.txt

    cat ${id}_MP_reads_Ig_Onco.txt ${id}_SA_reads_Ig_Onco.txt ${id}_MP_reads_Onco_Ig.txt ${id}_SA_reads_Onco_Ig.txt | \
    awk -v id=${id} -F'\t' '{print id,$0}' OFS='\t' > ${id}_tmp_reads.txt

    cut -f3 ${id}_tmp_reads.txt | sort | uniq > potential_reads.txt
    cat ${OncoBed} ${IgBed} > potential_regions.txt

    samtools view \
        -o "${id}_potential_reads.bam" \
        -F ${exclude_flag_split} \
        -M -L potential_regions.txt \
        -N potential_reads.txt \
        -X ${bam} ${bai}
    samtools index "${id}_potential_reads.bam"

    echo -e "ID\tREADTYPE\tQNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tRNEXT\tPNEXT\tTLEN\tSEQ\tQUAL" | cat - "${id}_tmp_reads.txt" > "${id}_mini_reads.txt"

    """
    outputs = {  "pseudo_sam" : "*_mini_reads.txt",
        "potential_reads_bam" : "*_potential_reads.bam",
        "potential_reads_bai" : "*_potential_reads.bam.bai"
              }
    docker = "jbalberge/samtools_cloud:1.13"


# ============================================================================
# STEP A: cluster_reads  (R container) -- clusters, finds breakpoints, ALT.
#   No REF counting here. Writes a breakpoints file for the samtools step.
# ============================================================================
class cluster_reads(wolf.Task):
    name = "cluster_reads"
    inputs = {
        "id":None,
        "OncoBed":None,
        "IgBed":None,
        "reads":None,
        "genome":"hg19",
        "eps":500,
        "minPts":2
    }
    script = """
    set -euxo pipefail

    Rscript - \
    "${reads}" \
    "${OncoBed}" \
    "${IgBed}" \
    "${genome}" \
    "${eps}" \
    "${minPts}" \
    "${id}_clusters_eps${eps}_minPts${minPts}_raw.txt" \
    "${id}_clusters_eps${eps}_minPts${minPts}_filtered.txt" \
    "${id}_clusters_eps${eps}_minPts${minPts}.pdf" \
    "${id}_breakpoints.txt" << "HEREDOC"

sem <- function(x) sd(x)/sqrt(length(x))

annotate.position <- function(hit.chr, hit.pos, bed) {
  name <- bed %>% filter(chr == hit.chr & start <= hit.pos & end >= hit.pos) %>% pull(name)
  if( identical(name, character(0))) {
    return("")
  } else {
    return(name)
  }
}

run.dbscan <- function(df, eps=1000, min.pts=2){
  mat <- df %>% select(POS.IG, POS.ONCO) %>% as.matrix()
  res <- dbscan(mat, eps, minPts = min.pts, weights = NULL, borderPoints = TRUE)
  df$Cluster=res$cluster
  df
}

# This script takes potential reads from CatchTheFISH,
# pairs reads, and cluster them, to call translocations

# INIT --------------------------------------------------------------------

set.seed(123)

library(tidyverse)
library(dbscan)
library(circlize)

args = commandArgs(trailingOnly=TRUE)

# INPUT
reads <- args[1]
onco.bed <- args[2]
ig.bed <- args[3]

# PARAMS
genome <- as.character(args[4])
eps <- as.numeric(args[5])
minPts <- as.numeric(args[6])

# OUTPUT
clustered.tx <- args[7]
clustered.tx.filtered <- args[8]
circos.plot <- args[9]
breakpoints.out <- args[10]

# breakpoints file always gets written (even empty) so the next task has input
bp.colnames <- c("ID","Chr_Cluster","SIDE","CHR","BP","PARTNER_CHR","PARTNER_BP","N_Molecules")
write.table(setNames(data.frame(matrix(ncol=length(bp.colnames), nrow=0)), bp.colnames),
            breakpoints.out, row.names = FALSE, quote = FALSE, sep = "\t")

# Parse input files -------------------------------------------------------

tx.header <- c("ID","SVTYPE","READNAME","FLAG","CHR1","POS1","QUAL","CIGAR","CHRB","POSB","TLEN","SEQ","SEQQUAL")
tx.col.classes <- c("character", "character", "character", "numeric", "character", "numeric", "numeric", "character", "character", "numeric", "numeric", "character", "character")
tx <- read.delim(reads, col.names = tx.header, colClasses = tx.col.classes, header = TRUE)

if( nrow(tx ) < 2 ) {
  res.colnames <- c("ID","Chr_Cluster","IG","ONCO","CHR.IG","IG_START","IG_END","CHR.ONCO","PARTNER_START","PARTNER_END","INTERVAL_LENGTH","IG_LENGTH","reads","maxMAPQ.ONCO","maxMAPQ.IG","N_Matches","N_Molecules","N_Mate_Pairs","N_Split_Reads","Breakpoint","BP.ONCO","BP.IG","BP_PRECISE","ALT")
  res <- setNames(data.frame(matrix(ncol = length(res.colnames), nrow = 0)), res.colnames)
  pdf(circos.plot, paper = "a4")
  circos.clear()
  circos.par("track.height" = 0.1, start.degree = 90)
  circos.initializeWithIdeogram(species = genome,  track.height = convert_height(2, "mm"), c("ideogram", "labels"))
  dev.off()
} else {

  bed.header <- c("chr", "start", "end", "name")
  bed.col.classes <- c("character", "numeric", "numeric", "character")
  onco.bed <- read.delim(onco.bed, col.names = bed.header, colClasses = bed.col.classes, header = FALSE)
  ig.bed <- read.delim(ig.bed, col.names = bed.header, colClasses = bed.col.classes, header = FALSE)

  tx.2 <- tx %>%
    rowwise() %>%
    mutate(IG=annotate.position(hit.chr = CHR1, hit.pos = POS1, bed = ig.bed),
           ONCO=annotate.position(hit.chr = CHR1, hit.pos = POS1, bed = onco.bed)) %>%
    group_by(ID, SVTYPE, READNAME) %>%
    filter(n()>=2)

  tx.ig <- tx.2 %>%
    filter(IG!="") %>%
    transmute(ID, SVTYPE, READNAME, CHR.IG=CHR1, POS.IG=POS1, IG, QUAL.IG=QUAL)

  tx.onco <- tx.2 %>%
    filter(ONCO!="") %>%
    transmute(ID, SVTYPE, READNAME, CHR.ONCO=CHR1, POS.ONCO=POS1, ONCO, QUAL.ONCO=QUAL)

  tx.join <- inner_join(tx.ig, tx.onco, by=c("ID", "SVTYPE", "READNAME")) %>%
    group_by(ID, SVTYPE, READNAME) %>%
    slice_max(n = 1, order_by = QUAL.ONCO, with_ties = FALSE)

  tx.join.annot <- tx.join %>%
    left_join(onco.bed, by=c("ONCO"="name"), suffix=c("", "oncoBed")) %>%
    left_join(ig.bed, by=c("IG"="name"), suffix=c("", "_igBed")) %>%
    mutate(trx_name=
             paste0("t(",
                    min(as.numeric(CHR.IG), as.numeric(CHR.ONCO)),
                    ";",
                    max(as.numeric(CHR.IG), as.numeric(CHR.ONCO)),
                    ")"),
           trx_full_name=paste0(IG, "-", ONCO))

  # Clustering --------------------------------------------------------------

  res1 <- tx.join.annot %>%
    group_by(ID, ONCO, IG) %>%
    group_modify(~ {
      distances <-  dist(matrix(c(.x$POS.IG, .x$POS.ONCO), ncol=2, byrow = FALSE), method = "maximum");
      .x %>% mutate(Cluster = dbscan(distances, eps=eps, minPts=minPts)$cluster)}) %>%
    mutate(Chr_Cluster=paste0(IG, "-", ONCO, "_", Cluster)) %>%
    filter(Cluster!=0) %>%
    group_by(ID, Chr_Cluster, IG, ONCO) %>%
    summarise(CHR.IG=first(CHR.IG),
              IG_START=min(POS.IG),
              IG_END=max(POS.IG),
              CHR.ONCO=first(CHR.ONCO),
              PARTNER_START=min(POS.ONCO),
              PARTNER_END=max(POS.ONCO),
              INTERVAL_LENGTH=PARTNER_END-PARTNER_START,
              IG_LENGTH=IG_END-IG_START,
              reads=paste(READNAME, collapse = ","),
              MAPQ.ONCO=paste(QUAL.ONCO,  collapse = ","),
              MAPQ.IG=paste(QUAL.IG,  collapse = ","),
              maxMAPQ.ONCO=max(QUAL.ONCO),
              maxMAPQ.IG=max(QUAL.IG),
              N_Matches=n(),
              N_Molecules=length(unique(READNAME)),
              N_Mate_Pairs=sum(SVTYPE=="Paired-Read"),
              N_Split_Reads=sum(SVTYPE=="Split-Read"),
              SA_ONCO=ifelse(any(SVTYPE=="Split-Read"),
                             median(POS.ONCO[SVTYPE=="Split-Read"]), NA_real_),
              SA_IG=ifelse(any(SVTYPE=="Split-Read"),
                           median(POS.IG[SVTYPE=="Split-Read"]), NA_real_),
              Breakpoint=paste0("Breakpoint_", first(Cluster))) %>%
    ungroup()

  # Breakpoint assignment: SA-precise if a split read exists, else span midpoint
  res1 <- res1 %>%
    mutate(
      BP_PRECISE = ifelse(N_Split_Reads > 0, "PRECISE", "IMPRECISE"),
      BP.ONCO = ifelse(N_Split_Reads > 0, round(SA_ONCO),
                       round((PARTNER_START + PARTNER_END)/2)),
      BP.IG   = ifelse(N_Split_Reads > 0, round(SA_IG),
                       round((IG_START + IG_END)/2)),
      ALT = N_Molecules
    )

  res <- res1 %>%
      filter(MAPQ.ONCO>=55 & INTERVAL_LENGTH>=1 & IG_LENGTH>=1)

  # Write breakpoints file (long format: one row per side per cluster) -------
  # Each row carries this side's CHR/BP AND the PARTNER side's CHR/BP, so the
  # downstream count_ref task can test whether a read's SA tag / discordant mate
  # points at the partner breakpoint (= ALT evidence for THIS junction).
  # N_Molecules is retained for comparison against the recounted ALT.
  if (nrow(res1) >= 1) {
    bp.onco <- res1 %>% transmute(ID, Chr_Cluster, SIDE="ONCO",
                                  CHR=CHR.ONCO, BP=BP.ONCO,
                                  PARTNER_CHR=CHR.IG, PARTNER_BP=BP.IG,
                                  N_Molecules=N_Molecules)
    bp.ig   <- res1 %>% transmute(ID, Chr_Cluster, SIDE="IG",
                                  CHR=CHR.IG, BP=BP.IG,
                                  PARTNER_CHR=CHR.ONCO, PARTNER_BP=BP.ONCO,
                                  N_Molecules=N_Molecules)
    bp.all  <- bind_rows(bp.onco, bp.ig) %>%
      filter(!is.na(BP) & BP >= 1 & !is.na(PARTNER_BP) & PARTNER_BP >= 1)
    write.table(bp.all, breakpoints.out, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  # Circos plot -------------------------------------------------------------

  pdf(circos.plot, paper = "a4")
  circos.clear()
  circos.par("track.height" = 0.1, start.degree = 90)
  circos.initializeWithIdeogram(species = genome,  track.height = convert_height(2, "mm"), c("ideogram", "labels"))

  if(nrow(res)>=1) {

    annot.partners <- res %>%
      transmute(chr=ifelse(startsWith(CHR.ONCO, "chr"), CHR.ONCO, paste0("chr", CHR.ONCO)),
                start=as.numeric(PARTNER_START),
                end=as.numeric(PARTNER_END),
                name=as.character(ONCO))

    annot.ig <- res %>%
      transmute(chr=ifelse(startsWith(CHR.IG, "chr"), CHR.IG, paste0("chr", CHR.IG)),
                start=as.numeric(IG_START),
                end=as.numeric(IG_END),
                name=as.character(IG))

    circos.genomicTrack(annot.partners, ylim=c(0, 0.1),
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, ytop = 0.1, ybottom = 0,...)
                          circos.genomicText(region, value, y=0, labels.column = 1, facing = "clockwise", niceFacing = TRUE, cex=0.8, adj=degree(0), ...)
                        }, bg.border = NA )

    circos.genomicLink(region1 = annot.partners,
                       region2 = annot.ig,
                       col = "red",
                       lwd = 10)
  }

  title(paste(unique(res$ID), collapse=", "))
  dev.off()
}

# Write output ------------------------------------------------------------

write.table(res, clustered.tx.filtered, row.names = FALSE, quote = TRUE, sep = "\t")
write.table(res1, clustered.tx, row.names = FALSE, quote = TRUE, sep = "\t")

HEREDOC

    """
    outputs = {
        "clusters_tx" : "*_raw.txt",
        "clusters_tx_filtered" : "*_filtered.txt",
        "circos_tx" : "*.pdf",
        "breakpoints" : "*_breakpoints.txt"
    }
    docker = "jbalberge/r-dbscan:latest"


# ============================================================================
# STEP B: count_ref  (samtools container) -- counts REF-spanning reads at each
#   breakpoint from the ORIGINAL bam, then computes VAF = ALT/(ALT+REF).
#   Runs on the same image as extract_sam, which has gcloud + samtools.
# ============================================================================
class count_ref(wolf.Task):
    name = "count_ref"
    inputs = {
        "id":None,
        "breakpoints":None,
        "gs_clean_bam":None,
        "gs_clean_bai":None,
        "bait_bed":None,
        "ref_min_mapq_onco":30,
        "ref_min_mapq_ig":0
    }
    overrides = {"gs_clean_bam":"string", "gs_clean_bai":"string"}
    script = """
    set -euxo pipefail

    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

    bam="${gs_clean_bam}"
    bai="${gs_clean_bai}"
    bp="${breakpoints}"
    baitbed="${bait_bed}"
    minq_onco="${ref_min_mapq_onco}"
    minq_ig="${ref_min_mapq_ig}"

    # ---- REF awk: concordant reads that strongly span the breakpoint ----
    # Drops SA-tagged (split) reads; counts reads whose CIGAR-implied span
    # strictly crosses BP. This is the REF (reference-allele) evidence.
    cat > count_span.awk << 'AWK'
function reflen_of(cig,  i,c,num,total){
  total=0; num="";
  for(i=1;i<=length(cig);i++){
    c=substr(cig,i,1);
    if(c ~ /[0-9]/){ num=num c }
    else { if(c=="M"||c=="D"||c=="N"||c=="="||c=="X") total+=num+0; num="" }
  }
  return total
}
BEGIN{ c=0 }
{
  if($0 ~ /\tSA:Z:/) next;
  start=$4; end=start+reflen_of($6)-1;
  if(start < BP && end > BP) c++
}
END{ print c }
AWK

    # Per-side table. ALT = N_Molecules (unique, detection-vetted molecule count).
    # REF = breakpoint-spanning reads, counted with a PER-SIDE mapq floor:
    #   ONCO side -> minq_onco (strict, unique-mapping locus)
    #   IG side   -> minq_ig   (loose, repetitive locus; matches detection mapq=0)
    # IN_BAIT = whether this breakpoint falls in a bait/probe interval.
    printf 'ID\tChr_Cluster\tSIDE\tCHR\tBP\tALT\tREF\tVAF\tIN_BAIT\n' > "${id}_vaf_byside.txt"

    # --- bait membership helper (chr-normalized, skips browser/track header) ---
    # Reads only the bait bed; query CHR/BP passed via -v; prints YES/NO in END.
    cat > bait_member.awk << 'AWK'
function norm(c){ sub(/^chr/,"",c); return c }
{
  if($1 ~ /^(browser|track)/) next;
  if($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) next;
  c=norm($1); n[c]++; s[c,n[c]]=$2+0; e[c,n[c]]=$3+0;
}
END{
  chr=norm(QCHR); pos=QBP+0; inb="NO";
  cnt=n[chr];
  for(i=1;i<=cnt;i++){ if(pos>=s[chr,i] && pos<=e[chr,i]){ inb="YES"; break } }
  print inb;
}
AWK

    # Breakpoints columns: ID Chr_Cluster SIDE CHR BP PARTNER_CHR PARTNER_BP N_Molecules
    tail -n +2 "${bp}" | while IFS=$'\t' read -r BID CLUST SIDE CHR BP PCHR PBP NMOL; do
        [ -z "${BP:-}" ] && continue

        region="${CHR}:${BP}-${BP}"

        # per-side mapq floor
        if [ "${SIDE}" = "ONCO" ]; then minq=${minq_onco}; else minq=${minq_ig}; fi

        # REF: clean spanning reads at this breakpoint, mapq >= minq
        REF=$(samtools view -f 2 -F 3852 -q ${minq} \
                -X ${bam} ${bai} "${region}" \
                | awk -v BP=${BP} -f count_span.awk)
        REF=${REF:-0}

        ALT=${NMOL:-0}

        # Per-side VAF = ALT / (ALT + REF). If denom is 0 -> NA.
        VAF=$(awk -v a="${ALT}" -v r="${REF}" 'BEGIN{ d=a+r; if(d>0) printf "%.6f", a/d; else print "NA" }')

        # bait membership for this breakpoint
        INB=$(awk -v QCHR="${CHR}" -v QBP="${BP}" -f bait_member.awk "${baitbed}")
        INB=${INB:-NO}

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${BID}" "${CLUST}" "${SIDE}" "${CHR}" "${BP}" "${ALT}" "${REF}" "${VAF}" "${INB}" >> "${id}_vaf_byside.txt"
    done

    # ---- Combine the two sides into one final VAF per cluster ----------------
    # Precedence:
    #   1. Bait membership (we trust high-coverage baited breakpoints):
    #        both sides in bait    -> mean(VAF.ONCO, VAF.IG)
    #        exactly one in bait   -> use that side's VAF
    #   2. Neither side in bait    -> fall back to order-of-magnitude rule:
    #        both real & same OOM  -> mean
    #        both real & diff OOM  -> bigger
    #   NA handling underlies all: a NA side can't be used; if the chosen side
    #   is NA, fall through to the other side / OOM as appropriate.
    cat > combine.awk << 'AWK'
function oom(x){ return (x>0) ? int(log(x)/log(10) + ( (log(x)/log(10)<0 && (log(x)/log(10))!=int(log(x)/log(10))) ? -1 : 0)) : "NA" }
BEGIN{ FS=OFS="\t" }
NR==1 { next }   # skip header
{
  cl=$2; side=$3; vaf=$8; inb=$9;
  clust_id[cl]=$1"\t"cl;
  if(side=="ONCO"){ onco[cl]=vaf; onco_b[cl]=inb }
  else if(side=="IG"){ ig[cl]=vaf; ig_b[cl]=inb }
}
END{
  print "ID","Chr_Cluster","VAF.ONCO","VAF.IG","ONCO_IN_BAIT","IG_IN_BAIT","VAF.FINAL","RULE";
  for(cl in clust_id){
    o=onco[cl]; g=ig[cl];
    ob=onco_b[cl]; gb=ig_b[cl];
    o_na=(o=="NA"||o==""); g_na=(g=="NA"||g=="");
    o_ok=(!o_na); g_ok=(!g_na);
    obait=(ob=="YES"); gbait=(gb=="YES");
    rule=""; final="NA";

    # 1. bait-based decision, but only trust a baited side if it's non-NA
    if(obait && gbait && o_ok && g_ok){ final=(o+g)/2.0; rule="both_bait_mean" }
    else if(obait && o_ok && !(gbait && g_ok)){ final=o; rule="onco_bait_only" }
    else if(gbait && g_ok && !(obait && o_ok)){ final=g; rule="ig_bait_only" }
    # 2. neither side usable via bait -> OOM fallback on whatever real values exist
    else if(o_ok && g_ok){
      oo=oom(o+0); og=oom(g+0);
      if(oo==og){ final=(o+g)/2.0; rule="fallback_same_oom_mean" }
      else { final=(o+0>=g+0)?o:g; rule="fallback_diff_oom_bigger" }
    }
    else if(o_ok){ final=o; rule="onco_only" }
    else if(g_ok){ final=g; rule="ig_only" }
    else { final="NA"; rule="both_NA" }

    of=(o_na?"NA":sprintf("%.6f",o+0));
    gf=(g_na?"NA":sprintf("%.6f",g+0));
    ff=(final=="NA"?"NA":sprintf("%.6f",final+0));
    print clust_id[cl], of, gf, ob, gb, ff, rule;
  }
}
AWK

    awk -f combine.awk "${id}_vaf_byside.txt" > "${id}_vaf.txt"
    """
    outputs = {
        "vaf" : "*_vaf.txt",
        "vaf_byside" : "*_vaf_byside.txt"
    }
    docker = "jbalberge/samtools_cloud:1.13"