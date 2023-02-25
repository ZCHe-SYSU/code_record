# Imports
import pandas as pd
import argparse
import os

# Own libraries
from data_manipulation.data import Data
from models.visualization.clusters import plot_cluster_images, plot_wsi_clusters, plot_wsi_clusters_interactions


##### Main #######
parser = argparse.ArgumentParser(description='Report cluster images from a given Leiden cluster configuration.')
parser.add_argument('--meta_folder',         dest='meta_folder',         type=str,            default=None,                   help='Purpose of the clustering, name of folder.')
parser.add_argument('--meta_field',          dest='meta_field',           type=str,            default=None,                  help='Meta field to use for the Logistic Regression or Cox event indicator.')
parser.add_argument('--matching_field',      dest='matching_field',       type=str,            default=None,                  help='Key used to match folds split and H5 representation file.')
parser.add_argument('--resolution',          dest='resolution',          type=float,           default=None,                  help='Minimum number of tiles per matching_field.')
parser.add_argument('--dpi',                 dest='dpi',                 type=int,            default=1000,                   help='Highest quality: 1000.')
parser.add_argument('--fold',                dest='fold',                type=int,            default=0,                      help='Minimum number of tiles per matching_field.')
parser.add_argument('--dataset',             dest='dataset',             type=str,            default='TCGAFFPE_LUADLUSC_5x', help='Dataset to use.')
parser.add_argument('--complete_path',       dest='complete_path',    type=str,            required=True,                     help='file path to run the leiden clustering folds.')
parser.add_argument('--additional_path',     dest='additional_path',  type=str,            default=None,                      help='Additional representation to assign leiden clusters.')
parser.add_argument('--min_tiles',           dest='min_tiles',           type=int,            default=400,                    help='Minimum number of tiles per matching_field.')
parser.add_argument('--dbs_path',            dest='dbs_path',            type=str,            default=None,                   help='Path for the output run.')
parser.add_argument('--img_size',            dest='img_size',            type=int,            default=224,                    help='Image size for the model.')
parser.add_argument('--img_ch',              dest='img_ch',              type=int,            default=3,                      help='Number of channels for the model.')
parser.add_argument('--marker',              dest='marker',              type=str,            default='he',                   help='Marker of dataset to use.')
parser.add_argument('--tile_img',            dest='tile_img',            action='store_true', default=False,                  help='Dump cluster tile images.')
parser.add_argument('--extensive',           dest='extensive',           action='store_true', default=False,                  help='Flag to dump test set cluster images in addition to train.')
parser.add_argument('--additional_as_fold',  dest='additional_fold',  action='store_true', default=False,                     help='Flag to specify if additional file will be used for cross-validation.')
args               = parser.parse_args()
meta_folder        = args.meta_folder
meta_field         = args.meta_field
matching_field     = args.matching_field
resolution         = args.resolution
dpi                = args.dpi
fold               = args.fold
min_tiles          = args.min_tiles
image_height       = args.img_size
image_width        = args.img_size
image_channels     = args.img_ch
marker             = args.marker
dataset            = args.dataset
complete_path      = args.complete_path
additional_path    = args.additional_path
dbs_path           = args.dbs_path
tile_img           = args.tile_img
extensive          = args.extensive
additional_fold    = args.additional_fold

# Dominating clusters to pull WSI.
only_id = True
value_cluster_ids = dict()
# value_cluster_ids[1] = []
# value_cluster_ids[0] = []


# Default path for GDC manifest.
manifest_csv = '%s/utilities/files/LUADLUSC/gdc_manifest.txt' % os.path.dirname(os.path.realpath(__file__))

# Default DBs path. 
if dbs_path is None:
    dbs_path = os.path.dirname(os.path.realpath(__file__))

# Leiden convention name.
groupby = 'leiden_%s' % resolution

# Dataset images.
data = Data(dataset=dataset, marker=marker, patch_h=image_height, patch_w=image_width, n_channels=image_channels, batch_size=64, project_path=dbs_path, load=True)

# Dump cluster images.
if tile_img:
    plot_cluster_images(groupby, meta_folder, data, fold, complete_path, dpi, value_cluster_ids, extensive=extensive)

# Save WSI overlay with clusters.
plot_wsi_clusters(groupby, meta_folder, matching_field, meta_field, data, fold, complete_path, additional_path, additional_fold, dpi, min_tiles, manifest_csv=manifest_csv,
                  value_cluster_ids=value_cluster_ids, type_='percent', only_id=only_id, n_wsi_samples=3)

# Save WSI overlay with clusters.
# plot_wsi_clusters_interactions(groupby, meta_folder, 'slides', meta_field, data, fold, complete_path, additional_path, additional, dpi, min_tiles, manifest_csv=manifest_csv,
#                                inter_dict=inter_dict, type_='percent', only_id=only_id, n_wsi_samples=2)