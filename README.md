# 3D-NaissIJ

3D-NAissIJ is an analysis image workflow to measure the presence of protein in the cytoplasm and membrane of cells in a fresh cardiac tissue sample prepared following the 3D-NaissI (3D-Native Tissue Imaging) method (https://doi.org/10.1007/s00018-025-05595-y)

3D-NaissIJ is implemented as a FIJI script and a Cellpose model. The main steps are :

- Cytoplasm segmentation : Cell cytoplasm are segmented using a dedicated cellpose model trained on WGA signal applied to each Z. Then segmentation are merged. They can be reduced in size by erosion and signal from other channel con be subtracted.
- The surface segmentation is derived from the cytoplasm segmentation by label dilation. Signal from other channel con be subtracted.
- Protein fraction is measured in cell surface and cytoplasm after thresholding with a specified threshold

## Getting started

### Installation

1. Install FIJI : https://fiji.sc/
2. Install the MorpholibJ plugin for FIJI : https://github.com/ijpb/MorphoLibJ
3. Install Cellpose v3.1 : https://github.com/MouseLand/cellpose
4. Install the PTBIOP plugin for FIJI : https://wiki-biop.epfl.ch/en/ipa/fiji/update-site
5. Download the 3D-NaissiJ cellpose model :
6. Download the 3D-NaissiJ FIJI script from this repositery

### Usage

1. Start FIJI and run the 3D-NaissiJ FIJI script with the script editor or by adding the script to the menu : https://imagej.net/scripting/
2. Choose script parameters :
   - Check save results to automatically save measurement and segemntation preview in the image folder
   - For the cytoplasme segmentation, corresponding parameters are :
     - CELLPOSE Python environment folder and model file : Specify the path for the Cellpose environment folder (containing python.exe) and the path for the 3D-NaissiJ cellpose model file.
     - CELLPOSE Channel and diameter : Specify the WGA channel as well as a cell diameter in pixel
     - MERGE LABELS, minimal IoU : Specify the a minimal Intersection Over Union to merge cellpose segmentation along the Z axis (around 0.6 - 0.8)
     - OPTIONNAL: labels erosion : Specify a value in pixel to reduce cytoplasm size
     - REMOVE SIGNAL Channel and Autothreshold method: Specify a channel which signal is removed from the cytoplasm. The autotheshold method will be applied to each Z on this channel. The advised autothreshold methods are Li, Moments or Otsu.
   - The surface segmentation is derived from the cytoplasm segmentation with the following parameters :
     - Label dilation : Specify the thickness of surface in pixel
     - REMOVE SIGNAL Channel and Autothreshold method: Specify a channel which signal is removed from the cytoplasm. The autotheshold method will be applied to each Z on this channel. The advised autothreshold methods are Li, Moments or Otsu.
   - Specify the Channels for the protein fraction measurements as well as a signal threshold

Note : The 3D-NaissiJ FIJI script is compatible with ImageJ Batch processing : https://imagej.net/scripting/batch

## Acknowledgments

If you use 3D-NaissIJ, please cite
