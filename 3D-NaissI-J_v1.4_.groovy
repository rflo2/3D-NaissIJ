// 3D-NaissI-J
// =================================

// Rémy Flores-Flores, INSERM

/*
 * Protein expression measurement in fresh cardiac tissue samples in the cytoplasme and cell surface
 * 
 * v1.1 07/03/2025
 * Corrected a resulttable bug when no surface, changed display image name, changed UI messages
 * 
 * v1.2 12/03/2025
 * added channel removal for surface, corrected bug : Li was always used for channel removal for cyto
 * 
 * v1.3 03/07/2025
 * work on files and batch mode
 * new segmentation display with surface border as overlay
 * cooccurence measurement on signal (manual thresholding, medin filter 1)
 * removed surface segmentation with dystrophin thresholding, surface segmentation with dilation from cytoplasm segmentation
 * 
 * v1.4 08/09/2025
 * Correction of surface segmentation bug
 * Measurement and signal segmentation code in functions
 * Add an optionnal second protein measurement aswell as co-occurence
 * Add autothreshold options for signal segmentation
 *  
 */

#@ File (label="Select a file") image_path
#@ Boolean (label="save results") do_save
#@ String (visibility=MESSAGE, value="CYTOPLASM SEGMENTATION", required=false) msgCS
#@ File (label="CELLPOSE: Python envionment folder",style="directory") python_folder
#@ File (label="CELLPOSE: model file",style="file") model_file
#@String(label="CELLPOSE: Channel", choices={"1","2","3","4","5","6","7","8","9"}) CP_ch
#@ Float (label="CELLPOSE: diameter",value=300) CP_dia
#@ Float (label="MERGE LABELS: minimal IoU",value=0.7) min_IoU
#@ Float (label="OPTIONNAL: labels erosion [pixel]",value=0) label_erosion
#@String(label="REMOVE SIGNAL: Channel (0: no effect)", choices={"0","1","2","3","4","5","6","7","8","9"}) channel_removal_cyto
#@ String (label="REMOVE SIGNAL: Autothreshold method",value="Li") method_cyto

#@ String (visibility=MESSAGE, value="SURFACE SEGMENTATION", required=false) msgDis
#@ Float (label="Label dilation [pixel]",value=20) label_dilation
#@String(label="REMOVE SIGNAL: Channel (0: no effect)", choices={"0","1","2","3","4","5","6","7","8","9"}) channel_removal_surf
#@ String (label="REMOVE SIGNAL: Autothreshold method",value="Li") method_surf

#@ String (visibility=MESSAGE, value="PROTEIN 1", required=false) msgME1
#@String(label="Channel", choices={"1","2","3","4","5","6","7","8","9"}) SI_ch1
#@ String (label="Fraction: Threshold",value=80) seuil1

#@ String (visibility=MESSAGE, value="PROTEIN 2, if 0 not used", required=false) msgME2
#@String(label="Channel", choices={"0","1","2","3","4","5","6","7","8","9"}) SI_ch2
#@ String (label="Fraction: Threshold",value=80) seuil2


//#@ String (visibility=MESSAGE, value="Segment protein, input value or autothreshold method : Li, Otsu, Moments", required=false) msgSD

import ij.IJ
import ij.ImagePlus
import ij.plugin.Duplicator
import ij.plugin.ImageCalculator
import ij.plugin.ZProjector
import ij.measure.ResultsTable
import ij.gui.TextRoi 
import ij.gui.Overlay
import loci.plugins.in.ImporterOptions
import loci.plugins.BF

import inra.ijpb.binary.distmap.ChamferMask2D
import inra.ijpb.label.filter.ChamferLabelErosion2DShort
import inra.ijpb.label.filter.ChamferLabelDilation2DShort
import inra.ijpb.binary.BinaryImages as BI
import inra.ijpb.label.LabelImages as LI
import inra.ijpb.math.ImageCalculator as IC
import inra.ijpb.measure.IntensityMeasures
import inra.ijpb.plugins.LabelToValuePlugin
import inra.ijpb.measure.region2d.Centroid

import java.awt.Font
import java.awt.Color
import org.apache.commons.io.FilenameUtils

IJ.run("Close All");

dir = image_path.getParent()
image_name = image_path.getName() 
image_basename = FilenameUtils.getBaseName( image_name )


imp_option = new ImporterOptions()
imp_option.setId(image_path.toString())
ImagePlus image = BF.openImagePlus(imp_option)[0]


wheat_germ = new Duplicator().run(image, CP_ch.toInteger(), CP_ch.toInteger(), 1, image.getNSlices(), 1, 1);

wheat_germ.show();
IJ.run("Cellpose ...", "env_path="+python_folder+" env_type=conda model= model_path=["+model_file+"] diameter="+CP_dia+" ch1=1 ch2=-1 additional_flags=--use_gpu");
cellpose = IJ.getImage();
cellpose.hide();
println("cellpose done")

cellpose_merged = mergeLabels(cellpose,min_IoU)
println("merge labels done")
cellpose_merged.show()
cellpose_eroded = erodeLabels(cellpose_merged,label_erosion)
println("erode labels done")
cellpose_eroded.show()


cellpose_dilated = dilateLabels(cellpose_eroded,label_dilation)
ImageCalculator.run(cellpose_dilated, cellpose_eroded, "Subtract stack");

if(channel_removal_cyto!="0"){
	removeSignalFromLabelsWithStackAutoThreshold(image,cellpose_eroded,channel_removal_cyto.toInteger(),method_cyto);
}
println("remove signals done")


if(channel_removal_surf!="0"){
	removeSignalFromLabelsWithStackAutoThreshold(image,cellpose_dilated,channel_removal_surf.toInteger(),method_surf);
}
println("remove signals done")

max_label = ZProjector.run(cellpose_eroded,"max").getProcessor().getStats().max;
cellpose_eroded.setDisplayRange(0,max_label)
cellpose_eroded.setTitle("cytoplasme")

max_label = ZProjector.run(cellpose_dilated,"max").getProcessor().getStats().max;
cellpose_dilated.setDisplayRange(0,max_label)
cellpose_dilated.setTitle("surface")


//SIGNAL SEGMENTATION
signal1 = new Duplicator().run(image, SI_ch1.toInteger(), SI_ch1.toInteger(), 1, image.getNSlices(), 1, 1);
signal1.setTitle("signal1")
segmentation1 = signal_segmentation(signal1,seuil1)
segmentation1.setDisplayRange(0,1)
segmentation1.setTitle("segmentation1")

if(SI_ch2!="0"){
	signal2 = new Duplicator().run(image, SI_ch2.toInteger(), SI_ch2.toInteger(), 1, image.getNSlices(), 1, 1);
	signal2.setTitle("signal2")
	segmentation2 = signal_segmentation(signal2,seuil2)
	segmentation2.setDisplayRange(0,1)
	segmentation2.setTitle("segmentation2")
	
	both = ImageCalculator.run(segmentation1, segmentation2, "Multiply create stack");
	both.setDisplayRange(0,1)
	both.setTitle("both")
}


//MEASUREMENT
(cyto_labels,cyto_volume) = getLabel_morpho(cellpose_eroded)
(surf_labels,surf_volume) = getLabel_morpho(cellpose_dilated)

p1_mean_cyto = measure(signal1,cellpose_eroded)
p1_mean_surf = measure(signal1,cellpose_dilated)
p1_frac_cyto = measure(segmentation1,cellpose_eroded)
p1_frac_surf = measure(segmentation1,cellpose_dilated)

if(SI_ch2!="0"){
	p2_mean_cyto = measure(signal2,cellpose_eroded)
	p2_mean_surf = measure(signal2,cellpose_dilated)
	p2_frac_cyto = measure(segmentation2,cellpose_eroded)
	p2_frac_surf = measure(segmentation2,cellpose_dilated)
	both_cyto = measure(both,cellpose_eroded)
	both_surf = measure(both,cellpose_dilated)
}


//SHOW RESULTS
rt = new ResultsTable(cyto_labels.size())
for(i=0;i<cyto_labels.size();i++){
	cellInd = cyto_labels[i]
	rt.setValue("Cellule",i,cellInd)
	j = surf_labels.findIndexOf{it==cellInd}
	rt.setValue("Cyto: volume",i,cyto_volume[i])
	if(j != -1){rt.setValue("Surf: volume",i,surf_volume[j])}
	rt.setValue("P1 mean cyto",i,p1_mean_cyto[i])
	if(j != -1){rt.setValue("P1 mean surf",i,p1_mean_surf[i])}
	rt.setValue("P1 frac cyto",i,p1_frac_cyto[i])
	if(j != -1){rt.setValue("P1 frac surf",i,p1_frac_surf[i])}
	if(SI_ch2!="0"){
		rt.setValue("P2 mean cyto",i,p2_mean_cyto[i])
		if(j != -1){rt.setValue("P2 mean surf",i,p2_mean_surf[i])}
		rt.setValue("P2 frac cyto",i,p1_frac_cyto[i])
		if(j != -1){rt.setValue("P2 frac surf",i,p2_frac_surf[i])}
		rt.setValue("both cyto",i,both_cyto[i])
		if(j != -1){rt.setValue("both surf",i,both_surf[i])}
	}
}
rt.show("results")
if(do_save){
	IJ.saveAs("Results", dir +"\\" + image_basename + "_results.csv");
}

//DISPLAY
display = draw2DBoundaries(wheat_germ,cellpose_dilated)
if(do_save){
	IJ.saveAs(display, "Tiff", dir +"\\" + image_basename + "_segmentation.tif");
}

if(SI_ch2!="0"){
	segmentation1.show()
	segmentation2.show()
	IJ.run(image, "Merge Channels...", "c1=segmentation1 c2=segmentation2 create");
	merge = IJ.getImage();
	segmentation1.close()
	segmentation2.close()
	display = draw2DBoundaries(merge,cellpose_dilated)
	if(do_save){
		IJ.saveAs(display, "Tiff", dir +"\\" + image_basename + "_coloc.tif");
	}
}



//FUNCTIONS

def measure(input,labels){
	input.show()
	labels.show()
	IJ.run("Intensity Measurements 2D/3D", "input="+input.getTitle()+" labels="+labels.getTitle()+" mean volume");
	rt = ResultsTable.getActiveTable()
	mean = rt.getColumn("Mean")
	input.hide()
	labels.hide()
	return mean
}

def getLabel_morpho(imp){
	imp.show()
	IJ.run("Analyze Regions 3D", "volume surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
	rt = ResultsTable.getActiveTable()
	volume = rt.getColumn("Volume")
	def labels = []
	nb = rt.getCounter()
	for(i=0;i<nb;i++){
		labels << rt.getLabel(i).toInteger()
	}
	imp.hide()
	return [labels,volume]
}


def signal_segmentation(signal,seuil){
	segmentation = signal.duplicate()
	try{
		seuil = seuil.toFloat();
		IJ.setRawThreshold(segmentation, seuil.toFloat(), 70000);
		IJ.run(segmentation, "Convert to Mask", "background=Dark black");
	}
	catch(NumberFormatException e){
		segmentation.show()
		IJ.run("Auto Threshold", "method="+seuil+" white stack");
		segmentation.hide()
	}
	IJ.run(segmentation, "Median...", "radius=1 stack");
	IJ.run(segmentation, "Divide...", "value=255 stack");
	return segmentation
}

def draw2DBoundaries(image,label_image){
	new_image = image.duplicate()
	new_image.show();
	for(s = 1; s < new_image.getNSlices()+1; s++){
		//draw boundaries
		new_image.setPosition(1,s,1)
		label_image.setPosition(s)
		slice = new ImagePlus("slice",label_image.getProcessor());
		IJ.run(slice, "Label Boundaries", "");
		bound = IJ.getImage()
		IJ.selectWindow(new_image.getTitle())
		IJ.run("Add Image...", "image="+bound.getTitle()+" x=0 y=0 opacity=100 zero");
		bound.close();
		//draw label number
		font = new Font("SansSerif", Font.PLAIN, 24);
		overlay = new Overlay()
		ip = label_image.getProcessor();
		labels = LI.findAllLabels(ip)
		centroids = Centroid.centroids(ip,labels)
		for(l = 0; l<labels.size();l++){
			//println(labels[l])//+centroids[l])
			roi = new TextRoi(centroids[l][0],centroids[l][1], labels[l].toString(), font);
			roi.setStrokeColor(new Color(1.00, 1.00, 1.00));
			overlay.add(roi)		
		}
		ips = new_image.getProcessor();
		ips.drawOverlay(overlay)
	/*label_stack = label_image.getImageStack()
	//stack = image.getImageStack()
	for(s = 1; s < stack.getSize()+1; s++){
		ip = stack.getProcessor(s);
		bound = LI.labelBoundaries(ip)
		overlay = new Overlay()
	}*/
	}
	return new_image
}


def createRGBandAddLabels(image){
	font = new Font("SansSerif", Font.PLAIN, 24);
	display = image.duplicate()
	IJ.run(display, "RGB Color", "");
	stack = image.getImageStack()
	stack_display = display.getImageStack()
	for(s = 1; s < stack.getSize()+1; s++){
		ip = stack.getProcessor(s);
		labels = LI.findAllLabels(ip)
		centroids = Centroid.centroids(ip,labels)
		overlay = new Overlay()
		bound = LI.labelBoundaries(ip)
		for(l = 0; l<labels.size();l++){
			//println(labels[l])//+centroids[l])
			roi = new TextRoi(centroids[l][0],centroids[l][1], labels[l].toString(), font);
			roi.setStrokeColor(new Color(1.00, 1.00, 1.00));
			overlay.add(roi)		
		}
		ips = stack_display.getProcessor(s);
		ips.drawOverlay(overlay)
		stack_display.setProcessor(ips,s);
	}	
	return display
}


def erodeLabels(labels,radius){
	new_labels = labels.duplicate()
	stack = new_labels.getImageStack()
	eroder = new ChamferLabelErosion2DShort(ChamferMask2D.BORGEFORS, radius)
	for(s = 1; s < stack.getSize()+1; s++){
		ip = stack.getProcessor(s);
		stack.setProcessor(eroder.process(ip),s)
	}
	return new_labels
}


def dilateLabels(labels,radius){
	new_labels = labels.duplicate()
	stack = new_labels.getImageStack()
	dilater = new ChamferLabelDilation2DShort(ChamferMask2D.BORGEFORS, radius)
	for(s = 1; s < stack.getSize()+1; s++){
		ip = stack.getProcessor(s);
		stack.setProcessor(dilater.process(ip),s)
	}
	return new_labels
}

def mergeLabels(labels,min_IoU){
	stack = labels.getImageStack().convertToFloat()
	max_label = stack.getProcessor(1).getStats().max
	for(s = 2; s < stack.getSize()+1; s++){
		//println(s)
		ip1 = stack.getProcessor(s-1);
		ip2 = stack.getProcessor(s);
		imp1 = new ImagePlus("imp1",ip1);
		imp2 = new ImagePlus("imp2",ip2);
		medians = new IntensityMeasures(imp1,imp2).getMedian()//.getColumn("Median")
		for(l = 0; l<medians.size();l++){
			IoU = measureIoU(ip1,medians.getValue("Median", l),ip2,l+1)
			//println((l+1)+" "+medians.getValue("Median", l)+" "+IoU)
			if(IoU<min_IoU){
				max_label=max_label+1
				medians.setValue("Median", l,max_label)
			}
		}
		new_ip = LabelToValuePlugin.process(imp2,medians,"Median").getProcessor()
		stack.setProcessor(new_ip,s)	
	}
	new_labels = new ImagePlus(labels.getTitle()+"_merged",stack)
	new_labels.setDisplayRange(0,max_label)
	IJ.run(new_labels, "Macro...", "code=[if (v != v) v = 0; ] stack");
	return new_labels
}

def removeSignalFromLabelsWithStackAutoThreshold(image,labels,channel,method){
	single_channel = new Duplicator().run(image, channel, channel, 1, image.getNSlices(), 1, 1);
	IJ.run(single_channel, "Auto Threshold", "method="+method+" white stack");
	IJ.run(single_channel, "Median...", "radius=1 stack");
	IJ.run(single_channel, "Invert", "stack");
	IJ.run(single_channel, "Divide...", "value=255 stack");
	return ImageCalculator.run(labels, single_channel, "Multiply stack");
}

def measureIoU(ip1,label1,ip2,label2){
	int[] lbl1 = [label1]
	int[] lbl2 = [label2]
	label1 = BI.binarize(LI.keepLabels(ip1, lbl1))
	label2 = BI.binarize(LI.keepLabels(ip2, lbl2))
	IC.Operation op_and = IC.Operation.AND;
	and = IC.combineImages(label1, label2, op_and);
	IC.Operation op_or = IC.Operation.OR;
	or = IC.combineImages(label1, label2, op_or);
	return and.getStats().mean / or.getStats().mean;
}