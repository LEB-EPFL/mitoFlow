/*Get initial information about the file and 
set the channel to the photoconverted one.
Needs a multistack image with a line ROI
For arch cell I use a line width of 100px*/
roiManager("reset");
date = "240125_"
name = getTitle();
nameNoExt = replace(name, ".tif", "");
Stack.getDimensions(width, height, channels, slices, frames);
roiManager("Add");

//Change this if your photoconverted channel is not 3!
run("Duplicate...", "duplicate channels=2");

//First iteration to initialise the csv file
roiManager("Select", 0);
Stack.setFrame(1);
run("Plot Profile");
Plot.showValues();
Table.renameColumn("X", "distance (µm)");
Table.renameColumn("Y", "t1");
selectWindow("Plot of " + nameNoExt + "-1");
close();
selectWindow(nameNoExt + "-1.tif");

/*Loop over remaining frames of the stack 
and append the y-data to the csv file.*/
for (i = 2; i <= frames; i++) {
	Stack.setFrame(i);
	run("Plot Profile");
	selectWindow("Plot of " + nameNoExt + "-1");
	Plot.getValues(x, y);
	selectWindow("Results");
	Table.setColumn("t" + i, y);
	selectWindow("Plot of " + nameNoExt + "-1");
	close();
}

//Save the file in the image directory as well as the ROI
saveAs("Results", File.directory + date + nameNoExt + "_lineProfileGreen" + ".csv");
run("Close");
selectWindow(nameNoExt + "-1.tif");
run("Close");

//second channel
selectWindow(name);
run("Duplicate...", "duplicate channels=3");

//First iteration to initialise the csv file
roiManager("Select", 0);
Stack.setFrame(1);
run("Plot Profile");
Plot.showValues();
Table.renameColumn("X", "distance (µm)");
Table.renameColumn("Y", "t1");
selectWindow("Plot of " + nameNoExt + "-1");
close();
selectWindow(nameNoExt + "-1.tif");

/*Loop over remaining frames of the stack 
and append the y-data to the csv file.*/
for (i = 2; i <= frames; i++) {
	Stack.setFrame(i);
	run("Plot Profile");
	selectWindow("Plot of " + nameNoExt + "-1");
	Plot.getValues(x, y);
	selectWindow("Results");
	Table.setColumn("t" + i, y);
	selectWindow("Plot of " + nameNoExt + "-1");
	close();
}

//Save the file in the image directory as well as the ROI
saveAs("Results", File.directory + date + nameNoExt + "_lineProfile" + ".csv");
run("Close");
roiManager("Save", File.directory + nameNoExt +"_lineProfile.roi");
close("*");