// input: directory containing subdirectories Stack 1, Stack 2, Stack 3, etc each containing:
//        -conventional image: .ome.tiff
//        -deconvolved image: _5ADVMLE.tif
//        -Masked Polycomb image: _MaskedPcBodies.tif
//        -Single cell maskes from conventional stack produced in TANGO: Nucleus 1.tif, Nucleus 2.tif etc
// NB: highly dependent on having the correct file structure. Stack folders must contain these and nothing else.
// output: Single cell PC images + 3D Objects Counter summary


//FUNCTION THAT CUTS SINGLE CELL POLYCOMB BODIES
///////////////////////////////////////////////
function action(inputImg) {
listImg = getFileList(inputImg);
for (i = 0; i < listImg.length; i++){              //for all Mask Pc images iamges in a folder...
if (endsWith(listImg[i], "_MaskedPcBodies.tif")){
					core = substring(listImg[i],0,lengthOf(listImg[i])-19);    
					open(inputImg + core + "_MaskedPcBodies.tif");
					PcImg=getTitle();
					
}}

for (n = 0; n < listImg.length; n++){            //for all conventional iamges in a folder...
			if (endsWith(listImg[n], ".ome.tiff")){
				core3 = substring(listImg[n],0,lengthOf(listImg[n])-9);
				open(inputImg + core3 + ".ome.tiff");
				Conv=getTitle();
}}		

for (n = 0; n < listImg.length; n++){            //for all single cell masks
			if (startsWith(listImg[n], "Nucleus")){	
				core2 = substring(listImg[n],0,lengthOf(listImg[n])-9);
				print(inputImg + "Nucleus " + n-2 +".tif");       // minus 2 as there are three other files (0,1,2) including conv stack, decon stack, bodies stack.
				open(inputImg + "Nucleus " + n-2 + ".tif");
				
				run("Invert", "stack");
				run("Multiply...", "value=200 stack");
				run("Multiply...", "value=200 stack");
				run("Multiply...", "value=200 stack");			// Invert and increase intensity values to ensure subtraction is effective
				SCMask=getTitle();
				//saveAs("Tiff", output+core2+"_PCBodies");
				
				
				
				imageCalculator("Subtract create stack", PcImg, SCMask);	// Subtract mask of nucleus to get z-stack with only Polycomb body signal of single nucleus
				Final=getTitle();
				saveAs("Tiff", inputImg+"Nucleus"+" "+ n-2 +"_SingleCellBodies");  // name it same as original

				imageCalculator("Subtract create stack", Conv, SCMask);	// Do the same but from the conventional file to get a z-stack with a whole single nucleus' signal
				Final2=getTitle();
				saveAs("Tiff", inputImg+"Nucleus"+" "+ n-2 +"_WholeNucleus");  // name it same as original - now have two files: one with whole nucleus, one with bodies, linked to the same number in the same stack
										
			}}			
       close();
     
}
///////////////////////////////////////////////



// FUNCTION THAT PERFORMS 3D OC AND SAVES .txt RESULTS FILES
////////////////////////////////////////////////////////////
function ObjectsCounter(inputImg) {
listImg = getFileList(inputImg);
for (i = 0; i < listImg.length; i++){              //for all Mask Pc images iamges in a folder...
if (endsWith(listImg[i], "_SingleCellBodies.tif")){
					core = substring(listImg[i],0,lengthOf(listImg[i])-21);    
					open(inputImg + core + "_SingleCellBodies.tif");
					SCImg=getTitle();
						run("3D OC Options", "volume integrated_density mean_gray_value median_gray_value dots_size=2 font_size=8 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
					    run("3D Objects Counter", "threshold=3 slice=21 min.=28 max.=20000000 exclude_objects_on_edges objects statistics summary");
						saveAs("txt",inputImg+core+"_bodies_summaryObjectsCounter");	// Run 3D Objects Counter with the given parameters. No upper limit is given so that objects can be filtered later. A text file is deposited in the stack folder corresponding to the same nucleus. This can be reviewed to check that everything is working satisfactorily.
						close(SCImg);	
					
}}

for (i = 0; i < listImg.length; i++){              //for all Mask Pc images iamges in a folder...
if (endsWith(listImg[i], "_WholeNucleus.tif")){
					core = substring(listImg[i],0,lengthOf(listImg[i])-17);    
					open(inputImg + core + "_WholeNucleus.tif");
					SCImg=getTitle();
						run("3D OC Options", "volume integrated_density mean_gray_value median_gray_value dots_size=2 font_size=8 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
					    run("3D Objects Counter", "threshold=3 slice=21 min.=28 max.=20000000 objects statistics summary");
						saveAs("txt",inputImg+core+"_nucleus_summaryObjectsCounter");	// Run 3D Objects Counter with the given parameters. This should identify a single object - the whole nucleus. A text file is deposited in the stack folder corresponding to the same nucleus.
						close(SCImg);	

}}
}
////////////////////////////////////////////////////////////


//----------------- LOOP THAT PERFORMS ACTION ON EVERY SUBDIRECTORY -------------------------------

input = getDirectory("Select directory:");	// Allow the user to select the appropriate directory.
list = getFileList(input);
print(list[0]);
setBatchMode(true);
for (i = 0; i < list.length; i++){   
if (endsWith(list[i], "/")){                //if is a subdirectory (Stack 1, Stack 2, ... ,)
					subDir = input + list[i]; 
                    action(subDir);
                    ObjectsCounter(subDir);
                    close("*");
}}
setBatchMode(false);	








