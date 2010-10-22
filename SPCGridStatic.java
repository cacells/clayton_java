/*
 * This is the main part of the simulation of the SPC version of the model 
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

import javax.swing.JOptionPane;
 
public class SPCGridStatic {

	public ArrayList <SPCCell> tissue;// List of cells that make up the tissue
	private Random rand = new Random();
	public static int maxlineage;//probably static is  unnecessary
	public int gsize;
	public static boolean writetogrid = false,precheck=true,migrateanywhere = true;;
	public static boolean smallcolony=false,forcinggrow=false;
	public int [] SCs,cellclonecount,clonesizecount;
	public int SCcount=0;
	public int countcolonies=0,countviable=0;
	double log2 = Math.log(2.0);

	
	public SPCGridStatic(int size, int maxC, double frac, boolean justSC) {//beth: this is the constructor
		// Create new instance of simulation with size of grid maximum SPC cycle and fraction of stem cells 
	    SPCCell.maxCycle = maxC+1;// (see SPCCell)
	    gsize = size;
		SPCBoxStatic[][] grid = new SPCBoxStatic[size][size];
        //beth: matrix of dimensions sizexsize containing homes, and called 'grid'
		// Temporary 2D array to hold boxes in Cartesian grid so that connections can be made

		tissue = new ArrayList<SPCCell>();// Creates the list structure for the cells that constitute the tissue
		if (smallcolony) makegroupofcells(grid,maxC,frac);
		else{
        if (justSC) placecolonies(grid,maxC,frac);
        else fillwithcells(grid,maxC,frac);
		}
        //note grid contents will be changed by method bcos I am passing a pointer to the grid 
        //although the grid pointer can't be changed (not that I want to!)
		//beth: at this point there are no holes (type 0 cells).
		for (int x = 0; x < size; x++) { //  Loop through all the boxes in the grid 
			for (int y = 0; y < size; y++) {
		        for (int xx = x - 1; xx <= x + 1; xx++) {
			        for (int yy = y - 1; yy <= y + 1; yy++) {
						if((y!=yy)||(x!=xx)) // Form links with their 8 immediate neighbours
			            grid[x][y].addNeighbour(grid[bounds(xx,size)][bounds(yy,size)]);
						//This maintains the cartesian relationship between each of the boxes without having to maintain the array
			        }
			    }
			}
	    } 
		//reset the values in the static cell counting arrays
		SPCCell.resetstaticcounters();
       if (writetogrid) {
    	   writeGrid();
    	   writetogrid = false;
       }
       SCs= new int[maxlineage];
       cellclonecount = new int[maxlineage];
       clonesizecount = new int[7];
	}//beth: end of constructor
	
	public SPCGridStatic(int size, int maxC, double frac, boolean justSC,String fname) {//constructor to read from file
		// Create new instance of simulation with size of grid maximum SPC cycle and fraction of stem cells 
	    SPCCell.maxCycle = maxC+1;// (see SPCCell)
	    gsize = size;
		SPCBoxStatic[][] grid = new SPCBoxStatic[size][size];
        //beth: matrix of dimensions sizexsize containing homes, and called 'grid'
		// Temporary 2D array to hold boxes in Cartesian grid so that connections can be made

		tissue = new ArrayList<SPCCell>();// Creates the list structure for the cells that constitute the tissue
        readGrid(grid,fname);
        //note grid contents will be changed by method bcos I am passing a pointer to the grid 
        //although the grid pointer can't be changed (not that I want to!)
		//beth: at this point there are no holes (type 0 cells).
		for (int x = 0; x < size; x++) { //  Loop through all the boxes in the grid 
			for (int y = 0; y < size; y++) {
		        for (int xx = x - 1; xx <= x + 1; xx++) {
			        for (int yy = y - 1; yy <= y + 1; yy++) {
						if((y!=yy)||(x!=xx)) // Form links with their 8 immediate neighbours
			            grid[x][y].addNeighbour(grid[bounds(xx,size)][bounds(yy,size)]);
						//This maintains the cartesian relationship between each of the boxes without having to maintain the array
			        }
			    }
			}
	    } 
		//reset the values in the static cell counting arrays
		SPCCell.resetstaticcounters();
		SCs= new int[maxlineage];
	       cellclonecount = new int[maxlineage];
	       clonesizecount = new int[7];
	}
	private void makegroupofcells(SPCBoxStatic[][] grid,int maxC,double frac){
		int x,y,k,radiussq = 10*10,cradiussq;
		int lineage =0;
		SPCCell cell;//just a name to use for each cell as it is placed in a home in the grid
		int sc,ncells=0;
		//make a cluster of PMB cells
		for (x = 0; x < gsize; x++) {
        	for(y=0;y<gsize;y++){
        		cradiussq = (x-31)*(x-31) + (y-31)*(y-31);
        		if (cradiussq < radiussq) {
        			ncells++;
			    grid[x][y] = new SPCBoxStatic(x,y);// New instance of SPCBox created and added to 2D grid
				cell = new SPCCell(grid[x][y],lineage);// New instance of SPCCell created and given unique lineage id
				grid[x][y].occupant = cell;// The new cell is added to the SPCBox
				cell.type=2;
				tissue.add(cell);// Add new cell to list of cells that constitute the tissue
				//lineage++;//uncomment to lineage all
        		}
        		else{
    			    grid[x][y] = new SPCBoxStatic(x,y);// New instance of SPCBox created and added to 2D grid
    				cell = new SPCCell(grid[x][y],lineage);// New instance of SPCCell created and given unique lineage id
    				grid[x][y].occupant = cell;// The new cell is added to the SPCBox
    				cell.type=0;//a space
    				tissue.add(cell);// Add new cell to list of cells that constitute the tissue
    				//lineage++;//uncomment to lineage all
        		}
			}
		}
		lineage++;//for only lineaging SC
		sc = (int)(ncells*frac);
		maxlineage = sc+1;//if only lineaging SC
		while(sc>0){ // while not enough stem cells allocated
			cell = tissue.get(rand.nextInt(tissue.size())); // Pick a cell at random from tissue list
			if (cell.type==2){ // If that cell is not already a stem cell type 1  
				cell.type=1;// Chance to an SC
				cell.lineage = lineage;
				lineage++;
				sc--;
			}
		}

		
	}
	private void fillwithcells(SPCBoxStatic[][] grid,int maxC,double frac){
		int x,y,k;
		int lineage =0;
		SPCCell cell;//just a name to use for each cell as it is placed in a home in the grid
		int sc = (int)(gsize*gsize*frac);
		for (k = 0; k < gsize; k++) {
			x = k;
        	for(y=0;y<k;y++){
			    grid[x][y] = new SPCBoxStatic(x,y);// New instance of SPCBox created and added to 2D grid
				cell = new SPCCell(grid[x][y],lineage);// New instance of SPCCell created and given unique lineage id
				grid[x][y].occupant = cell;// The new cell is added to the SPCBox
				cell.type=2;
				tissue.add(cell);// Add new cell to list of cells that constitute the tissue
				//lineage++;//uncomment to lineage all
			    grid[y][x] = new SPCBoxStatic(y,x);// New instance of SPCBox created and added to 2D grid
				cell = new SPCCell(grid[y][x],lineage);// New instance of SPCCell created and given unique lineage id
				grid[y][x].occupant = cell;// The new cell is added to the SPCBox
				cell.type=2;
				tissue.add(cell);// Add new cell to list of cells that constitute the tissue
				//lineage++;//uncomment to lineage all
			}
		    grid[k][k] = new SPCBoxStatic(k,k);// New instance of SPCBox created and added to 2D grid
			cell = new SPCCell(grid[k][k],lineage);// New instance of SPCCell created and given unique lineage id
			grid[k][k].occupant = cell;// The new cell is added to the SPCBox
			cell.type=2;
			tissue.add(cell);// Add new cell to list of cells that constitute the tissue
			//lineage++;//uncomment to lineage all
		}
		lineage++;//for only lineaging SC
		maxlineage = sc+1;//if only lineaging SC
		while(sc>0){ // while not enough stem cells allocated
			cell = tissue.get(rand.nextInt(tissue.size())); // Pick a cell at random from tissue list
			if (cell.type!=1){ // If that cell is not already a stem cell type 1  
				cell.type=1;// Chance to an SC
				cell.lineage = lineage;
				lineage++;
				sc--;
			}
		}

		
	}
	
	private void placecolonies(SPCBoxStatic[][] grid,int maxC,double frac){
		int x,y,k;
		int lineage = 0;
		SPCCell cell;//just a name to use for each cell as it is placed in a home in the grid
		int sc = (int)(gsize*gsize*frac);//beth: 64 hard coded.could use size
		//System.out.println("sc "+sc);
		for (k = 0; k < gsize; k++) {
			x = k;
        	for(y=0;y<k;y++){
			    grid[x][y] = new SPCBoxStatic(x,y);// New instance of SPCBox created and added to 2D grid
				cell = new SPCCell(grid[x][y],lineage);// New instance of SPCCell created and given unique lineage id
				grid[x][y].occupant = cell;// The new cell is added to the SPCBox
				cell.type=0;// space
				tissue.add(cell);// Add new cell to list of cells that constitute the tissue
			    grid[y][x] = new SPCBoxStatic(y,x);// New instance of SPCBox created and added to 2D grid
				cell = new SPCCell(grid[y][x],lineage);// New instance of SPCCell created and given unique lineage id
				grid[y][x].occupant = cell;// The new cell is added to the SPCBox
				cell.type=0;// Cell type set randomly to either SPC_1,SPC_2,SPC_3 or SPC_4
				tissue.add(cell);// Add new cell to list of cells that constitute the tissue
			}
		    grid[k][k] = new SPCBoxStatic(k,k);// New instance of SPCBox created and added to 2D grid
			cell = new SPCCell(grid[k][k],lineage);// New instance of SPCCell created and given unique lineage id
			grid[k][k].occupant = cell;// The new cell is added to the SPCBox
			cell.type=0;// Cell type set randomly to either SPC_1,SPC_2,SPC_3 or SPC_4
			tissue.add(cell);// Add new cell to list of cells that constitute the tissue
		}
		maxlineage = sc+1;//stem cells and spaces
		while(sc>0){ // while not enough stem cells allocated
			cell = tissue.get(rand.nextInt(tissue.size())); // Pick a cell at random from tissue list
			if (cell.type!=1){ // If that cell is not already a stem cell type 1 
				lineage++;
				cell.type=1;// Chance to an SC
				cell.lineage = lineage;
				//System.out.println("lineage "+lineage);
				sc--;
			}
		}

	}

	private int bounds(int a,int size) {  // Creates the toroidal links between top and bottom and left and right
		if (a < 0) return size + a;
		if (a >= size) return a - size;
		return a;
	}
	
	public void stain(){ // Stains all cells in the tissue list
		for (SPCCell c : tissue) {
		    if(c.type>0){
				c.stain=1.0;
			}
	    }
	}

	public void iterate() { // The main iterative loop of the simulation
		//beth: 

        SPCCell cHold;
		// Create a list to hold cells that are spaces or have the capacity to detach
        ArrayList<SPCCell> growArray = new ArrayList<SPCCell>();
        for (SPCCell c : tissue) { // loop through the tissue (ArrayList of cells)
		    c.maintain(); // Calls each cell to maintain its state re: detach and/or grow
			if(c.type==1)growArray.add(c); // If cell is an SC, add to grow list
			//if(c.canDetach)growArray.add(c);// If cell can detach add to grow list
	    }
        //beth: go through the list and see if anything grows into those spots
		while(growArray.size()>0){ // Randomly loop through the grow list
			cHold=growArray.remove(rand.nextInt(growArray.size()));
			cHold.grow();// Test to see if cell can be replaced by new proliferation
		}
	}
	public void iterateandcount(boolean printit) {
		//beth: 
		int[] colonies = new int[maxlineage];
		int logval,csc;
		boolean needstogrow;

        SPCCell cHold;
		// Create a list to hold cells that are spaces or have the capacity to detach
        ArrayList<SPCCell> growArray = new ArrayList<SPCCell>();
        
		SCcount=0;
		countcolonies=0;
		countviable=0;
        for (int ii=0;ii<maxlineage;ii++) {
        	SCs[ii] = 0;
        	cellclonecount[ii] = 0;
        }
        for (int ii=0;ii<7;ii++)clonesizecount[ii] = 0;
        
        for (SPCCell c : tissue) { // loop through the tissue (ArrayList of cells)
			colonies[c.lineage] = colonies[c.lineage] | 1;
			cellclonecount[c.lineage]++;
			if (c.type==1) {
				growArray.add(c);
				SCs[c.lineage] = SCs[c.lineage] | 1;
				SCcount++;
				//System.out.println("yo");
			}     	
		    c.maintainandcount(); // Calls each cell to maintain its state re: detach and/or grow
			//if(c.type==1)growArray.add(c); // If cell is an SC, add to grow list
			//if(c.canDetach)growArray.add(c);// If cell can detach add to grow list
	    }
		for(int ii=0;ii<maxlineage;ii++){
			countcolonies = countcolonies+ colonies[ii];
			countviable = countviable+ SCs[ii];
			csc = cellclonecount[ii];
			if (csc > 1){
				logval = (int) (Math.log(csc-1)/log2);
				if (logval > 6) logval = 6;
                clonesizecount[logval]++;
			}
		}
		if (printit) System.out.print(" "+countviable+"  "+countcolonies+"  "+SCcount);
        //beth: go through the list and see if anything grows into those spots
		while(growArray.size()>0){ // Randomly loop through the grow list
			cHold=growArray.remove(rand.nextInt(growArray.size()));
			if (precheck)
				cHold.growandcountprecheck(forcinggrow);
			else
				cHold.growandcount();// Test to see if cell can be replaced by new proliferation

		}
		if (printit) {
			System.out.print("  "+SPCCell.totalproliferations);

		for (int ii=0;ii<7;ii++){
			System.out.print(" "+clonesizecount[ii]);
		}
		System.out.println();
		}
	}

	public void iterateandcount_random(boolean printit) {
		//beth: 
		int[] colonies = new int[maxlineage];
		int logval,csc;
		boolean needstogrow;

        SPCCell cHold;
		// Create a list to hold cells that are spaces or have the capacity to detach
        ArrayList<SPCCell> growArray = new ArrayList<SPCCell>();
        ArrayList<SPCCell> migrateArray = new ArrayList<SPCCell>();
        
		SCcount=0;
		countcolonies=0;
		countviable=0;
        for (int ii=0;ii<maxlineage;ii++) {
        	SCs[ii] = 0;
        	cellclonecount[ii] = 0;
        }
        for (int ii=0;ii<7;ii++)clonesizecount[ii] = 0;
        
        for (SPCCell c : tissue) { // loop through the tissue (ArrayList of cells)
			colonies[c.lineage] = colonies[c.lineage] | 1;
			cellclonecount[c.lineage]++;
			if (c.type==1) {
				growArray.add(c);
				SCs[c.lineage] = SCs[c.lineage] | 1;
				SCcount++;
				//System.out.println("yo");
			}
			else {
				migrateArray.add(c);
			}
		    c.maintainandcount(); // Calls each cell to maintain its state re: detach and/or grow
			//if(c.type==1)growArray.add(c); // If cell is an SC, add to grow list
			//if(c.canDetach)growArray.add(c);// If cell can detach add to grow list
	    }
		for(int ii=0;ii<maxlineage;ii++){
			countcolonies = countcolonies+ colonies[ii];
			countviable = countviable+ SCs[ii];
			csc = cellclonecount[ii];
			if (csc > 1){
				logval = (int) (Math.log(csc-1)/log2);
				if (logval > 6) logval = 6;
                clonesizecount[logval]++;
			}
		}
		if (printit) System.out.print(" "+countviable+"  "+countcolonies+"  "+SCcount);
        //beth: go through the list and see if anything grows into those spots
		while(growArray.size()>0){ // Randomly loop through the grow list
			cHold=growArray.remove(rand.nextInt(growArray.size()));
			if (migrateanywhere) 
				cHold.growandcountrandom(migrateArray);
			else {
			if (precheck)
				cHold.growandcountprecheck(forcinggrow);
			else
				cHold.growandcount();// Test to see if cell can be replaced by new proliferation
			}

		}
		if (printit) {
			System.out.print("  "+SPCCell.totalproliferations);

		for (int ii=0;ii<7;ii++){
			System.out.print(" "+clonesizecount[ii]);
		}
		System.out.println();
		}
	}

public void writeGrid(){
	try{
		BufferedWriter bufGrid = new BufferedWriter(new FileWriter("grid.dat"));
		
        SPCCell cHold;
        for (SPCCell c : tissue) { // loop through the tissue (ArrayList of cells)
			bufGrid.write(c.home.x+" "+c.home.y+" "+c.type+" "+c.lineage);
			bufGrid.newLine();	
	    }

		bufGrid.close();
		System.out.println("Finished writing grid");
	}
	catch(IOException e){
	}
}


public void readGrid(SPCBoxStatic[][] grid, String fname){
	      try{
		  		// Open the file
		  		FileInputStream fstream = new FileInputStream(fname);
		  		// Get the object of DataInputStream
		  		DataInputStream in = new DataInputStream(fstream);
		          BufferedReader br = new BufferedReader(new InputStreamReader(in));
		  		String strLine;
		  		//Read File Line By Line
		  		int i1,i2,i3,i4;
				int fi = 0,li;
				int len;
				int count = 0;
		        SPCCell cell;
                maxlineage = 0;
                
		  		while ((strLine = br.readLine()) != null) 	{
		  			// Print the content on the console
		  			fi = 0;
		  			//System.out.println (strLine);
		  			len = strLine.length();
						li = strLine.indexOf(' ');
						i1 = Integer.parseInt(strLine.substring(fi,li));
						fi = li+1;
						li = strLine.indexOf(' ',fi);
						i2 = Integer.parseInt(strLine.substring(fi,li));
						fi = li+1;
						li = strLine.indexOf(' ',fi);
						i3 = Integer.parseInt(strLine.substring(fi,li));
						fi = li+1;
						i4 = Integer.parseInt(strLine.substring(fi,len));
						grid[i1][i2] = new SPCBoxStatic(i1,i2);
						cell = new SPCCell(grid[i1][i2],i4);
						if (maxlineage < i4) maxlineage = i4;
						grid[i1][i2].occupant = cell;
						cell.type = i3;
						tissue.add(cell);// Add new cell to list of cells that constitute the tissue
		  		}
		  		maxlineage++;
		  		//System.out.println("max lin "+maxlineage);
		  		//Close the input stream
		  		in.close();
		  		}catch (Exception e){//Catch exception if any
		  			System.err.println("Error: " + e.getMessage());
		  		}

}

}//end SPCGridStatic