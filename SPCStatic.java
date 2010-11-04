/*
 * This is a class that contains main and
 * will call the SPCGrid simulation and graphically display results.
 * usage: SPCStatic scfrac scrate maxiters smallcolony forcegrow readgrid
 * typical usage: SPCStatic 0.22 0.2 100 small blank new
 * for small group of cells, 0.22 frac SC, 5 iters per SC prolif, 100 iters total, 
 * don't force the grow, and generate a new grid
 * SPC is for Single Progenitor Cell
 */


import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.*;

import java.util.*;
 
public class SPCStatic extends JFrame implements Runnable {

    SPCGridStatic experiment;
    Random rand = new Random();    
	Thread runner;
    Container mainWindow;
	CAImagePanel CApicture;
	Image backImg1;
	Graphics backGr1;
	Image backImg2;
	Graphics backGr2;
	JProgressBar progressBar;
	JTextArea msgBtn,msgBtn2,msgBtn3;
	JPanel buttonHolderhigh;
	int scale = 20;//beth: could set to 1. Makes the colour transitions better?
	int border = 20;
	int iterations;
	int gSize;
    Colour palette = new Colour();
	int[] colorindices = {45,1,54};//{0,1,2,54,4,5};
	int nnw = colorindices.length-1;
//    Color[] colours = {Color.white,Color.black,Color.green,Color.blue,Color.yellow,Color.red,Color.pink};
    Color[] javaColours;
    double[][] epsColours;
    //Color[] colours = {Color.black,Color.white,Color.green,Color.blue,Color.yellow,Color.red,Color.pink};
    boolean writeImages = false,writetoGrid = true;
    static int maxIters = 300,writefreq=10;
	int lin=64*64;//max number of different cell lines equals cells in grid
	double dlin = (double) lin;
	public static boolean readgrid = false;
	BufferedWriter bufGrid;


	public SPCStatic(int size, int maxC, double frac) {
	    gSize=size;
	    if (readgrid)
			experiment = new SPCGridStatic(size, maxC, frac,false,"grid.dat");
	    else
		    experiment = new SPCGridStatic(size, maxC, frac,false);

		lin = SPCGridStatic.maxlineage;//in case the grid is not full
		dlin = (double) lin;
		

		setVisible(true);
		backImg1 = createImage(scale * size, scale * size);
		backGr1 = backImg1.getGraphics();
		backImg2 = createImage(scale * size, scale * size);
		backGr2 = backImg2.getGraphics();
		setpalette();
		
	    int wscale = 6;//scale for main panel
	    int btnHeight = 480-384;//found by trial and error - must be a better way!
	    //although no buttons yet
	    int wh = (gSize*1)*wscale + 2*border;// +btnHeight;//mainWindow height
	    int ww = (gSize*2)*wscale + 3*border;//mainWindow width   
	    
		mainWindow = getContentPane();
		mainWindow.setLayout(new BorderLayout());
		setSize(ww,wh);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setVisible(true);
        CApicture = new CAImagePanel(ww,wh);
        CApicture.setBorder(border);
        CApicture.rowstoShow = gSize;
        mainWindow.add(CApicture,BorderLayout.CENTER);
		setVisible(true);
		
		buttonHolderhigh = new JPanel();
		buttonHolderhigh.setLayout(new GridLayout(1,4));
		
        msgBtn = new JTextArea("   Colony count: "+0);
        msgBtn.setEditable(false);
        buttonHolderhigh.add(msgBtn);
        msgBtn2 = new JTextArea("   Viable Colony count: "+0);
        msgBtn2.setEditable(false);
        buttonHolderhigh.add(msgBtn2);
        msgBtn3 = new JTextArea("   Stem Cell count: "+0);
        msgBtn3.setEditable(false);
        buttonHolderhigh.add(msgBtn3);
		
	    progressBar = new JProgressBar(JProgressBar.HORIZONTAL,0,maxIters-1);
	    progressBar.setValue(0);
	    progressBar.setStringPainted(true);
	    
	    buttonHolderhigh.add(progressBar);	    
	    mainWindow.add(buttonHolderhigh, BorderLayout.SOUTH);
		setVisible(true);
		

		
	}
	//new ones
	public void drawCA() {

		int a;
		for (SPCCell c : experiment.tissue){
			a = c.type;
			CApicture.drawCircleAt(c.home.x, c.home.y, javaColours[a], 1);
		}
	    CApicture.updateGraphic();
	}
	public void drawCAstain() {
		double cstain;
  	    CApicture.clearCAPanel(2);
		for (SPCCell c : experiment.tissue){
			//if ((c.stain < minstain) && (c.type>0)) minstain = c.stain;
			cstain = c.stain;
			if (c.type==1) 
				CApicture.drawCircleAt2(c.home.x, c.home.y, palette.Javashades(cstain), 2);
			else 
				CApicture.drawCircleAt(c.home.x, c.home.y, palette.Javashades(cstain), 2);
		}
	    //outputImage();
	    CApicture.updateGraphic();
	}

	public void drawCALineage() {
		double celllin;
  	    CApicture.clearCAPanel(2);
  	    Color Lineagecolour;

		for (SPCCell c : experiment.tissue){

			//if ((c.stain < minstain) && (c.type>0)) minstain = c.stain;
			if (c.type==0){
				CApicture.drawCircleAt(c.home.x, c.home.y, javaColours[0], 2);
			}
			else{
				if (experiment.SCs[c.lineage]==1){
					celllin = (double)c.lineage/(dlin+1.0);
					if (c.type==1) {
						CApicture.drawCircleAt2(c.home.x, c.home.y, palette.Javagrey(celllin), 2);
					}
					else 
						CApicture.drawCircleAt(c.home.x, c.home.y, palette.Javagrey(celllin), 2);
				}
				else{//draw it blue - not viable
					celllin = (double)c.lineage/(dlin+1.0);
					if (c.type==1) {
						CApicture.drawCircleAt2(c.home.x, c.home.y, javaColours[2], 2);
					}
					else 
						CApicture.drawCircleAt(c.home.x, c.home.y, javaColours[2], 2);
				}
			}
		}

		//System.out.println("v "+countviable+" c "+countcolonies+" s "+SCcount+" p "+SPCCell.totalproliferations);
	    //outputImage();
	    CApicture.updateGraphic();
	}
	
	public void openoutputfile(){
		try{
			bufGrid = new BufferedWriter(new FileWriter("allgrid.dat"));
		}
		catch(IOException e){
			}
	}
	public void closeoutputfile(){
		try{
			bufGrid.close();
			System.out.println("Finished writing grid");
		}
		catch(IOException e){
			}
	}
	public void writeGrid(){
		try{

	        for (SPCCell c : experiment.tissue) { // loop through the tissue (ArrayList of cells)
				bufGrid.write(c.home.x+" "+c.home.y+" "+c.type+" "+c.age+" "+c.lineage);
				bufGrid.newLine();	
		    }


		}
		catch(IOException e){
		}
	}

	public void initialise(){
			CApicture.setScale(gSize,gSize,scale,gSize,gSize,scale);
      	    CApicture.clearCAPanel(1);
      	    CApicture.clearCAPanel(2);
      	    CApicture.clearParent();
		    iterations=0;
	}
	
	
	public void start() {
		initialise();
		if (runner == null) {
			runner = new Thread(this);
		}
		runner.start();
	}


	public void run() {
		if (runner == Thread.currentThread()) {
			drawCA();
			//either draw colonies:
			drawCALineage();
			System.out.println("scRate "+SPCCell.scRate);
			if (writeImages) CApicture.writeImage(0);
			openoutputfile();
			for(iterations=0; iterations<maxIters; iterations++){
				if(iterations==0)experiment.stain();// stain all cells at start

				experiment.iterateandcount_random(false);
				progressBar.setValue(iterations);
				if (iterations%writefreq==0)writeGrid();
				if (iterations%2==0){
					drawCA();
					//either draw colonies:
					drawCALineage();
					//or draw stain:
					//drawCAstain();
					msgBtn.setText("   Colony count: "+ experiment.countcolonies);
					msgBtn2.setText("   Viable Colony count: "+ experiment.countviable);
					msgBtn3.setText("   Stem Cell count: "+experiment.SCcount);
					if (writeImages) CApicture.writeImage(iterations+2);
				}
				SPCCell.resetstaticcounters();
				//if((iterations%5)==0)postscriptPrint("SPC"+iterations+".eps");
				// This will produce a postscript output of the tissue
			}
			closeoutputfile();
		}
	}


	public void postscriptPrint(String fileName) {
		int xx;
		int yy;
		int state;
		boolean flag;
		try {
			java.io.FileWriter file = new java.io.FileWriter(fileName);
			java.io.BufferedWriter buffer = new java.io.BufferedWriter(file);
			System.out.println(fileName);
			buffer.write("%!PS-Adobe-2.0 EPSF-2.0");
			buffer.newLine();
			buffer.write("%%Title: test.eps");
			buffer.newLine();
			buffer.write("%%Creator: gnuplot 4.2 patchlevel 4");
			buffer.newLine();
			buffer.write("%%CreationDate: Thu Jun  4 14:16:00 2009");
			buffer.newLine();
			buffer.write("%%DocumentFonts: (atend)");
			buffer.newLine();
			buffer.write("%%BoundingBox: 0 0 300 300");
			buffer.newLine();
			buffer.write("%%EndComments");
			buffer.newLine();
			for (SPCCell c : experiment.tissue){
				if(c.type>0){
					xx = (c.home.x * 4) + 20;
					yy = (c.home.y * 4) + 20;
					if (c.proliferated) {
						buffer.write("newpath " + xx + " " + yy + " 1.5 0 360 arc fill\n");
						buffer.write("0 setgray\n");
						buffer.write("newpath " + xx + " " + yy + " 1.5 0 360 arc  stroke\n");
					} else {
						buffer.write("0.75 setgray\n");
						buffer.write("newpath " + xx + " " + yy + " 1.5 0 360 arc fill\n");
					}
				}
			}
			buffer.write("showpage");
			buffer.newLine();
			buffer.write("%%Trailer");
			buffer.newLine();
			buffer.write("%%DocumentFonts: Helvetica");
			buffer.newLine();
			buffer.close();
		} catch (java.io.IOException e) {
			System.out.println(e.toString());
		}
	}
    public void setpalette(){
    	int ind = colorindices.length;
    	javaColours = new Color[ind];
    	epsColours = new double[ind][3];
    	for (int i=0;i<ind;i++){
    		//System.out.println("color index "+colorindices[i]);
    		javaColours[i] = palette.chooseJavaColour(colorindices[i]);
    		epsColours[i] = palette.chooseEPSColour(colorindices[i]);
    	}
    }

	public static void main(String args[]) {
		double initalSeed = 0.1;
		int arglen = args.length;
		maxIters = 1000;
		SPCStatic s;
		switch(arglen){
		case 7:
			if (args[6].contains("read")) readgrid = true;
		case 6:
			if (args[5].contains("random")) SPCGridStatic.migrateanywhere = true;
		case 5:
			if (args[4].contains("force")) SPCGridStatic.forcinggrow = true;
		case 4:
			if (args[3].contains("small")) {
				SPCGridStatic.smallcolony = true;
				readgrid = false;
			}
		case 3:
			maxIters = Integer.parseInt(args[2]);
		case 2:
			SPCCell.scRate = Double.parseDouble(args[1]);
            if (SPCCell.scRate > 1.0) SPCCell.scRate = 1.0;
		case 1:
			initalSeed = Double.parseDouble(args[0]);
			s = new SPCStatic(64, 1, initalSeed);
			break;
		default:
			writefreq = (int) Math.round(1.0/SPCCell.scRate);
			System.out.println("writefreq "+writefreq);
			maxIters = 100*writefreq;
			//maxiters set above, scrate set to 1 in spccell,force and small set in gridstatic.
			SPCGridStatic.migrateanywhere = true;
			SPCGridStatic.forcinggrow = true;
			s = new SPCStatic(64, 1, 0.22);
		}

		s.start();
//		if(arglen>0){
//			initalSeed = Double.parseDouble(args[0]);
//			if (arglen > 1) {
//				SPCCell.scRate = Double.parseDouble(args[1]);
//                if (SPCCell.scRate > 1.0) SPCCell.scRate = 1.0;
//			}
//			else SPCCell.scRate = 1.0;
//			if (arglen > 2) {
//				maxIters = Integer.parseInt(args[2]);
//			}
//			else maxIters = 100;
//			//System.out.println("scRate "+SPCCell.scRate);
//			SPCStatic s = new SPCStatic(64, 1, initalSeed);
//			s.start();
//		}else{
//			SPCStatic s = new SPCStatic(64, 1, 0.22);
//			s.start();
//		}
		
	}
}

