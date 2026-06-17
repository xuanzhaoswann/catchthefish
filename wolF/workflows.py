from .tasks import *

import wolf
from wolf.fc import SyncToWorkspace

def catchthefish(id,
                 gs_clean_bam,
                 gs_clean_bai,
                 OncoBed="gs://fc-15c5f539-a3ac-4ccd-8a9c-66bdd72cb5ff/catchthefish/ref/hg19_oncogene_regions.bed",
                 IgBed="gs://fc-15c5f539-a3ac-4ccd-8a9c-66bdd72cb5ff/catchthefish/ref/hg19_ig_regions.bed",
                 min_map_quality_onco=30,
                 min_map_quality_ig=0,
                 exclude_flag_split=1540,
                 exclude_flag_mate=1548,
                 include_flag_mate=1,
                 genome="hg19",
                 eps=500,
                 minPts=2,
                 add_gs_suffix=False,
                 workspace=None):

    bam_in_the_cloud = extract_sam(inputs={
                                        "id":id,
                                        "gs_clean_bam":gs_clean_bam,
                                        "gs_clean_bai":gs_clean_bai,
                                        "OncoBed":OncoBed,
                                        "IgBed":IgBed,
                                        "min_map_quality_onco":min_map_quality_onco,
                                        "min_map_quality_ig":min_map_quality_ig,
                                        "exclude_flag_split":exclude_flag_split,
                                        "exclude_flag_mate":exclude_flag_mate,
                                        "include_flag_mate":include_flag_mate,
                                        "add_gs_suffix":add_gs_suffix
                                      })

    tx = cluster_reads(inputs = {
        "id":id,
        "OncoBed":OncoBed,
        "IgBed":IgBed,
        "reads":bam_in_the_cloud["pseudo_sam"],
        "genome":genome,
        "eps":eps,
        "minPts":minPts
    })

    if workspace:
        upload_dict = {
                'clusters_tx':tx["clusters_tx"],
                'clusters_tx_filtered':tx["clusters_tx_filtered"],
                'clusters_candidate_bam':bam_in_the_cloud["potential_reads_bam"],
                'clusters_candidate_bai':bam_in_the_cloud["potential_reads_bai"],
                }
        sync_run = SyncToWorkspace(
                nameworkspace = workspace,
                #entity_type = "pair",
                entity_type = "sample",
                entity_name = id,
                attr_map = upload_dict
                )

