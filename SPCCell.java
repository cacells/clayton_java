/*
 * This contains the classes for the Cell object for the SPC version of the model
 */

import java.util.*;
 
class SPCCell{
	public static int maxCycle = 2; // number of different cell types
	public static Random rand = new Random();
	public static double rval = .08;
	public static double TAFraction = 1.0 - 2.0*rval;// The probability of a producing EP and PMB
	static Integer[] neighbours = {0,1,2,3,4,5,6,7};

	public boolean canDetach, canGrow,proliferated,isavail;
	public int type; // 0 = space, 1 = EP, 2=PMB
	public SPCBoxStatic home;// The box the cell sits in
	public double stain;
	public int age;
	public static double scRate=0.33;// Relative SC proliferation rate if scRate = 0.5 SC proliferation rate would be half SPC rate
	public int lineage;
	
	public static int[] cellcounts=new int[maxCycle+1];
	public static int[] prolifcounts=new int[maxCycle+1];
	public static double[] stainsums=new double[maxCycle+1];
	public static double[] agesums=new double[maxCycle+1];
	public static int totalproliferations=0;
	
	public SPCCell(SPCBoxStatic home,int lin){
		this.home=home;
		lineage = lin;
		canDetach=false;
		canGrow=false;
		proliferated=false;
		isavail = false;
		stain = 0.0;
		age = 0;
	}
	
	public static void resetstaticcounters(){
		for (int i=0;i<maxCycle+1;i++){
			cellcounts[i]=0;
			prolifcounts[i]=0;
			stainsums[i]=0.0;
			agesums[i]=0.0;
			totalproliferations = 0;
		}
	}
	
	public void maintain(){// Determines if a Cell can detach or grow 
		canDetach=(type==maxCycle);// For standard SPC model only PMB can detach
		canGrow = (type==1);// For standard SPC model only EP can grow
		isavail = (canDetach || (type==0));//PMB or space
		proliferated = false;
	}
	
	public void maintainandcount(){// Determines if a Cell can detach or grow and counts
		cellcounts[type]++;
		stainsums[type] = stainsums[type] + stain;
		agesums[type] = agesums[type] + age;
		canDetach=(type==maxCycle);// For standard SPC model only PMB can detach
		canGrow = (type==1);// For standard SPC model only EP can grow
		isavail = (canDetach || (type==0));//PMB or space
		proliferated = false;
	}
	
	public void growth(SPCCell cHold){ // Growth occurs into cell. chold is the neighbour that is taking over
		if(cHold.type==1){// only SC can proliferate
			if(rand.nextDouble()<TAFraction){ // The progenitor cell will be a PMB and parent stays as EP
				type=2;
			}else{ // if not
				if(rand.nextDouble()<0.5){ // then there is an equal probability that the cell will be
					type=1; //another EP
				}else{
					type=2;  // or that both cells will be PMB
					cHold.type=2;
				}
			}
		}		
			canGrow=false;// New cell will not proliferate again in this iteration
			cHold.canGrow=false;// Proliferating cell will not proliferate again in this iteration
			cHold.proliferated=true;// The proliferating cell has proliferated
			cHold.stain = cHold.stain/2.0;// Divide the label resting cell between the two cells
			stain = cHold.stain;// As above
			cHold.age = cHold.age + 1;// Divide the label resting cell between the two cells
			age = cHold.age;// As above
			lineage = cHold.lineage;// New cell takes on lineage of proliferating cell
			isavail = false;
	}
	public void growthandcount(SPCCell cHold){
		totalproliferations++;
		prolifcounts[cHold.type]++;
		if(cHold.type==1){// only SC can proliferate
			if(rand.nextDouble()<TAFraction){ // The progenitor cell will be a PMB and parent stays as EP
				type=2;
			}else{ // if not
				if(rand.nextDouble()<0.5){ // then there is an equal probability that the cell will be
					type=1; //another EP
				}else{
					type=2;  // or that both cells will be PMB
					cHold.type=2;
				}
			}
		}
		
		canGrow=false;// New cell will not proliferate again in this iteration
		cHold.canGrow=false;// Proliferating cell will not proliferate again in this iteration
		cHold.proliferated=true;// The proliferating cell has proliferated
		cHold.stain = cHold.stain/2.0;// Divide the label resting cell between the two cells
		stain = cHold.stain;// As above
		cHold.age = cHold.age + 1;// Divide the label resting cell between the two cells
		age = cHold.age;// As above
		lineage = cHold.lineage;// New cell takes on lineage of proliferating cell 
		isavail = false;
	}
	

	public boolean grow(){//new version of grow
		int sizeA = 8;//always 8 neighbours - otherwise this method needs changing
		ArrayList<Integer> nlist = new ArrayList<Integer>(Arrays.asList(neighbours));//initialise nlist
        int a,b;
		SPCCell cHold;
        for(int i=0;i<sizeA;i++){ // Loop from starting point through list of neighbours
        	a = rand.nextInt(nlist.size());//pick random list index
        	b = nlist.remove(a);//use the value at that index and make the list smaller
			cHold = home.getNeighbour(b);
			if(cHold.isavail){ // If neighbour is available for takeover
				if (rand.nextDouble()<scRate){//does we prolif?
					cHold.growth(this);// proliferate
					return true;// and stop search
				}
				else {
					canGrow = false;//don't prolif this round
					return false;
				}
			}
        }
		return false;// Return false if no proliferating cell can be found 
	}
	public boolean growandcount(){
		int sizeA = 8;//always 8 neighbours - otherwise this method needs changing
		ArrayList<Integer> nlist = new ArrayList<Integer>(Arrays.asList(neighbours));//initialise nlist
        int a,b;
		SPCCell cHold;
        for(int i=0;i<sizeA;i++){ // Loop from starting point through list of neighbours
        	a = rand.nextInt(nlist.size());//pick random list index
        	b = nlist.remove(a);//use the value at that index and make the list smaller
			cHold = home.getNeighbour(b);
			if(cHold.isavail){ // If neighbour is available for takeover
				//debug System.out.print("1 ");
				if (rand.nextDouble()<scRate){//do we prolif?
					//debug System.out.println("1 ");
					cHold.growthandcount(this);// proliferate
					return true;// and stop search
				}
				else {
					//debug System.out.println("0 ");
					canGrow = false;//don't prolif this round
					return false;
				}
			}
        }
		return false;// Return false if no proliferating cell can be found 
	}
	public boolean growandcountprecheck(boolean forcinggrow){
		int sizeA = 8;//always 8 neighbours - otherwise this method needs changing
		ArrayList<Integer> nlist = new ArrayList<Integer>(Arrays.asList(neighbours));//initialise nlist
		int a,b;
		SPCCell cHold;
		boolean found = true;
		if (rand.nextDouble()<scRate){//do we prolif?
			//find an available neighbour or grow upwards
			for(int i=0;i<sizeA;i++){ // Loop from starting point through list of neighbours
				a = rand.nextInt(nlist.size());//pick random list index
				b = nlist.remove(a);//use the value at that index and make the list smaller
				cHold = home.getNeighbour(b);

				if(cHold.isavail){ // If neighbour is available for takeover
					cHold.growthandcount(this);// proliferate
                    return true;
				}
			}
			//does it need to grow still, but there was no space?
			//if (forcinggrow) growup();
			if (forcinggrow){
				found = searchfurther(1);
			}
			return found;
		}	
		return false;// Return false if no proliferating cell can be found 
	}
	
	public boolean growandcountrandom(ArrayList <SPCCell> migrateArray){
		SPCCell cHold;
		boolean found = false;
		if (rand.nextDouble()<scRate){//do we prolif?
			//find an available neighbour or grow upwards
			while (!found && migrateArray.size()>0){
				cHold = migrateArray.remove(rand.nextInt(migrateArray.size()));
				if(cHold.isavail){ // If neighbour is available for takeover
					cHold.growthandcount(this);// proliferate
					found = true;
					return found;
				}
			}
		}
		return found;
	}
	
	private void growup(){//grow upwards because there is no space here...
		totalproliferations++;
		prolifcounts[type]++;
		//if(type==1){// only SC can proliferate
		/*uncomment for normal prolif rule.
		 * if (rand.nextDouble()<=rval) {
			type=2;//changes to a PMB
			isavail = false;//it probably was this already
		}*/
		canGrow=false;// Proliferating cell will not proliferate again in this iteration
		proliferated=true;// The proliferating cell has proliferated
		stain = stain/2.0;// Divide the label resting cell between the two cells
		age = age + 1;// Divide the label resting cell between the two cells
	}
	private boolean searchfurther(int n){
		//System.out.print("forcing ");
		//find a spare PMB amongst the neighbours of your neighbours
		ArrayList<Integer> nlist = new ArrayList<Integer>(Arrays.asList(neighbours));//initialise nlist
		ArrayList<Integer> nofnlist = new ArrayList<Integer>(Arrays.asList(neighbours));
		int sizeA = 8;//always 8 neighbours - otherwise this method needs changing
		int a,b,a2,b2;
		SPCCell cHold,cHold2;
		for(int i=0;i<sizeA;i++){ // Loop from starting point through list of neighbours
			a = rand.nextInt(nlist.size());//pick random list index
			b = nlist.remove(a);//use the value at that index and make the list smaller
			cHold = home.getNeighbour(b);
			nofnlist = new ArrayList<Integer>(Arrays.asList(neighbours));
			for(int ii=0;ii<sizeA;ii++){ // Loop from starting point through list of neighbours			
				a2 = rand.nextInt(nofnlist.size());//pick random list index
				b2 = nofnlist.remove(a2);//use the value at that index and make the list smaller
				cHold2 = cHold.home.getNeighbour(b2);		
                if(cHold2.isavail){ // If neighbour is available for takeover
				    cHold2.growthandcount(this);// proliferate
				    //System.out.println("... found");
                    return true;
                }
			}
		}
		return false;
	}
	private int bounds(int a,int size) {  // Creates the toroidal links between top and bottom and left and right
		if (a < 0) return size + a;
		if (a >= size) return a - size;
		return a;
	}
}
