filename = getTitle();
run("Make Inverse");
getDimensions(width, height, channels, slices, frames);

for (i = 1; i<=slices; i++){
	Stack.setSlice(i);
	for (j=1; j<=channels; j++){
		Stack.setChannel(j);
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
	}
}

save("Z:/_Lab members/Julius/03_DataAnalysis/MitoFlow/FISH_analysis/data_bg_cleared/" + filename);