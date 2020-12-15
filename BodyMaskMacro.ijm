// input: directory containing subdirectories Stack 1, Stack 2, Stack 3 etc. each containing:
//        -conventional image .ome.tiff
//        -deconvolved image _5ADVMLE.tif
// output: Masked polycomb bodies from conventional image saved as stack in Stack 1, Stack 2, ..., subdirectories
// NB: this script is highly dependent on having the prescribed file structure. If any additional files are present, it may not run correctly.

// Defined function to produce Polycomb body masks:

function action(inputImg) {
listImg = getFileList(inputImg);
for (i = 0; i < listImg.length; i++){              // for all deconvolved images in a folder...
if (endsWith(listImg[i], "_5ADVMLE.tif")){
					core = substring(listImg[i],0,lengthOf(listImg[i])-12);
					open(inputImg + core + "_5ADVMLE.tif");	// open the deconvolved file
					Decon=getTitle();

					// Make mask:
						run("Subtract Background...", "rolling=4 stack");
						setAutoThreshold("Otsu dark stack");
						run("Convert to Mask", "method=Otsu background=Dark");
						run("Invert", "stack");
						run("16-bit");
						run("Multiply...", "value=300 stack"); // Increases non-body region intensity to ensure subtraction removes all other signal
						Decon=getTitle();
}}

			for (n = 0; n < listImg.length; n++){            // for all conventional iamges in a folder...
			if (endsWith(listImg[n], ".ome.tiff")){
				core2 = substring(listImg[n],0,lengthOf(listImg[n])-9);
				open(inputImg + core2 + ".ome.tiff");
				Conv=getTitle();	// Just need title to save mask file with
			}}			
						
// Perform masking:
imageCalculator("Subtract create stack", Conv, Decon);	// Remove non-body regions from whole stack
Final=getTitle();
saveAs("Tiff", inputImg+core2+"_MaskedPcBodies");	// Save masked bodies file
       close();
}
///////////////////////////////////////////////////////////////////////////////////////////

// LOOP THAT PERFORMS ACTION ON EVERY SUBDIRECTORY

input = getDirectory("Select desired directory:"); // Allow the user to select the appropriate directory.
list = getFileList(input);
setBatchMode(true);
for (i = 0; i < list.length; i++){   
if (endsWith(list[i], "/")){                //if is a subdirectory (Stack 1, Stack 2, ... ,)
					subDir = input + list[i]; 
                    action(subDir);
                    close();
}}
setBatchMode(false);					



