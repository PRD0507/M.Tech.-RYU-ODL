import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;


public class IntervalCluster {

	private static final int NO_PARENT = -1; 
	int clusterNumber = 3;
	int TOTAL_DATA;
	int r = 2;
	double minX, minY;
	double maxX, maxY;
	double xLineDistance;
	double regionInterval;
	int pathCounter = 0;
	double shortestPathInsideCluster = 0;
	double bandwidth;
	
	double totalDistance;
	double insideClusterDistance;
	double outsideClusterDistance;
	
	List<Record> data = new ArrayList<Record>();
	List<Cluster> clusters = new ArrayList<Cluster>();
	List<ClusterHeadSet> clusterHeadList = new ArrayList<ClusterHeadSet>();
	double[][] adjacencyMatrix;
	
	Scanner scan = new Scanner(System.in);

	Map<Integer, List<Record>> clusterRecords = new HashMap<Integer, List<Record>>();
	static double execution;
	
	
	static int dataSize=100;
	String filename = "100x.txt";
	String filename2 = "100y.txt";
	double[] x= new double[dataSize];
	double[] y= new double[dataSize];
	
	public static void main(String args[]) {
		
		double sum=0;
		int sim=1;
	
		double[] executionTime = new double[sim];
		for(int k=0;k<sim;k++) {
		long start = System.nanoTime();
		
		IntervalCluster demo = new IntervalCluster();
		demo.readCSV();
		demo.genereateRecord();
		demo.initiateClusterAndCentroid();
		demo.printClusterInformation();
		
		demo.clusterHead();
		demo.printClusterHeadInformation();
		demo.createClusterHeadAdjaceny();
		demo.printClusterHeadAdjaceny();
		demo.shortestPathAlgo();
		
		long end = System.nanoTime();

		execution = end - start;
		System.out.println("Execution time: " + execution + " nanoseconds");
		
		double elapsedTimeInSecond = (double) execution / 1_000_000_000;
		System.out.println("Execution time: " + elapsedTimeInSecond + " seconds");

		elapsedTimeInSecond = elapsedTimeInSecond*Math.pow(10, 3);
		System.out.println("Execution time: " +  elapsedTimeInSecond + " milliseconds");
		
		
		executionTime[k]=elapsedTimeInSecond;
		sum=sum+elapsedTimeInSecond;
		demo.printResult();	
				
		}
		
		System.out.print("[");
		for(int k=0;k<sim;k++) {
			System.out.print(executionTime[k]+",");
		}
		System.out.print("]");
		
		System.out.println(" ");
		double avg=sum/sim;
		System.out.print("Average"+avg);
		//demo.result();
		
		
		}
	
	public void readCSV() {
		
		  try {
		      File myObj = new File("C:\\Users\\prdpr\\Desktop\\DataPoints\\"+filename);
		      Scanner myReader = new Scanner(myObj);
		      while (myReader.hasNextLine()) {
		        String data = myReader.nextLine();
		      //  System.out.println(data);
		        
		        String [] str = data.split("\\s+");
		        for(int i=0; i<str.length;i++) {
		        	 //System.out.println("i"+i+"str"+str[i]);
		        	 x[i]= Double.parseDouble(str[i]);
		        	 
		        	 
		        }
		        
		        
		      }
		      myReader.close();
		    } catch (FileNotFoundException e) {
		      System.out.println("An error occurred.");
		      e.printStackTrace();
		    }
		  
		  try {
		      File myObj = new File("C:\\Users\\prdpr\\Desktop\\DataPoints\\"+filename2);
		      Scanner myReader = new Scanner(myObj);
		      while (myReader.hasNextLine()) {
		        String data = myReader.nextLine();
		        //System.out.println(data);
		        
		        String [] str = data.split("\\s+");
		        for(int i=0; i<str.length;i++) {
		        	 //System.out.println("i"+i+"str"+str[i]);
		        	 y[i]= Double.parseDouble(str[i]);
		        	 
		        	 
		        }
		        
		        
		      }
		      myReader.close();
		    } catch (FileNotFoundException e) {
		      System.out.println("An error occurred.");
		      e.printStackTrace();
		    }
		  
		  
	}
	
	private void genereateRecord() {
		//First three record mark the centroids
		//Record record = new Record(1.0, 1.0,0);
		/*data.add(record);
		record = new Record(2.5, 3.5,0);
		data.add(record);
		record = new Record(5.0, 7.0,0);
		data.add(record);*/
		
		//Record Data
		/*Record record = new Record(0.0, 1.0,1);
		data.add(record);
		record = new Record(1.5, 2.0,2);
		data.add(record);
		record = new Record(1.25, 1.6,3);
		data.add(record);
		record = new Record(3.0, 4.0,4);
		data.add(record);
		record = new Record(5.0, 7.0,5);
		data.add(record);
		record = new Record(3.5, 5.0,6);
		data.add(record);
		record = new Record(4.5, 5.0,7);
		data.add(record);
		record = new Record(3.5, 4.5,8);
		data.add(record);
		record = new Record(6.5, 6,9);
		data.add(record);
		record = new Record(9.0, 9.0,10);
		data.add(record);
		record = new Record(7.0, 10.0,11);
		data.add(record);*/
		
		int j=1;
		for(int i=0;i<dataSize;i++) {
			
			Record record = new Record(x[i], y[i], j);
			data.add(record);
			j++;
		}
		
	}
	
	private void initiateClusterAndCentroid() {
			
		Iterator<Record> iterator = data.iterator();
		//Iterator<Centroid> citerator = cData.iterator();
		TOTAL_DATA = data.size();
		System.out.println("TOTAL_DATA"+TOTAL_DATA);
		
		
		System.out.println("Enter the min clustering interval");
		//regionInterval = scan.nextDouble();
		regionInterval = 30;
		
		
		Record rec = null;
		rec = iterator.next();
		minX = rec.getX();
		maxX = rec.getX();
		
		minY = rec.getY();
		maxY = rec.getY();
		
		while(iterator.hasNext()) {
			
			rec = iterator.next();
			
			if(rec.getX() <= minX) {
				minX = rec.getX();
			}
			if(rec.getX() >= maxX) {
				maxX = rec.getX();
			}
			if(rec.getY() <= minY)
			{
				minY = rec.getY();
			}
			if(rec.getY() >= maxY){
				maxY = rec.getY();
			}
		}
	//	System.out.println("min"+minX+"max"+maxX+"miny"+minY+"maxy"+maxY);
		xLineDistance = Math.sqrt(Math.pow(minX-maxX, 2) + Math.pow(0,2));
		clusterNumber = (int) Math.ceil(xLineDistance/regionInterval);
		
		System.out.println("cluster number"+clusterNumber);
		
		adjacencyMatrix = new double[clusterNumber][clusterNumber] ;
		
		System.out.println(xLineDistance);
		//System.out.println(regionInterval);
		
		int len = data.size();
		System.out.println("len"+len);
		int cNumDetected[] = new int[clusterNumber+1];
		
		for(int i=0; i<=clusterNumber;i++) {
			cNumDetected[i]=0;
		}
		
		for (int i = 0; i < len; i++) {
			
			//System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%");
			//System.out.println("I"+i);
			Record record = data.get(i);
			int cNumber = (int) Math.ceil(record.getX()/regionInterval);
			
		//	System.out.println("Cnu "+cNumber);
		
			if(cNumDetected[cNumber] == 0) {
			//	System.out.println("First Time ");
				
				if( cNumber == 0) {
				record.setClusterNumber(1);
				initializeCluster(1, record);
				cNumDetected[cNumber]++;
				cNumDetected[cNumber+1]++;
				}
				else {
				record.setClusterNumber(cNumber);
				initializeCluster(cNumber, record);
				cNumDetected[cNumber]++;
				}
			}
			else {
				
			//	System.out.println("Region Detected Before ");
				
				if( cNumber == 0) {
					record.setClusterNumber(1);
					clusterRecords.get(1).add(record);
				}
				else {
				
					record.setClusterNumber(cNumber);		
					clusterRecords.get(cNumber).add(record);
				}
				
				/*System.out.println("** Cluster Information **");
				for(Cluster cluster : clusters) {
					System.out.println(cluster);
				}
				System.out.println("*********************");*/
			}
		}	
	}
	
	private void initializeCluster(int cNumber, Record record) {
		
		Cluster cluster = new Cluster(cNumber);
		clusters.add(cluster);
		List<Record> clusterRecord = new ArrayList<Record>();
		clusterRecord.add(record);
		
		//System.out.println("Inside initializeCluster");
		//System.out.println("CNUM"+cNumber);
		clusterRecords.put(cNumber, clusterRecord);

	}
	
	private void printClusterInformation() {
		   System.out.println("****** FINAL CLUSTER INFORMATIN *********");
		   for (Map.Entry<Integer, List<Record>> entry : clusterRecords.entrySet())  {
	        System.out.println("Key = " + entry.getKey() + 
	                         ", Value = " + entry.getValue()); 
		   }
	}
	
	
	private void printClusterHeadInformation() {
		   System.out.println("****** Each Cluster Head INFORMATION *********");
		   for(ClusterHeadSet clusterHead : clusterHeadList) {
			   System.out.println(clusterHead);
		   }
	   }

	private void clusterHead() {
	
		double regionXCenter;
		regionXCenter = minX+(regionInterval/2);
		 
		for (Map.Entry<Integer, List<Record>> entry : clusterRecords.entrySet())  {
			
			int counter=0;
			double minDistClusterArray[] = new double[TOTAL_DATA];
			double x[] = new double[TOTAL_DATA];
			double y[] = new double[TOTAL_DATA];
			double nodeNumber[] = new double[TOTAL_DATA];
			
	        Integer key = entry.getKey();
	        List<Record> recordList = entry.getValue();
	        
	        
	        int length = recordList.size();
	        
	        double max = recordList.get(0).getY();
	        double min = recordList.get(0).getY();
	        
	        for(int i=1; i<length;i++) {
	        	
	        	if(max < recordList.get(i).getY()) {
	        		max = recordList.get(i).getY();
	        	}  	
	        }
	        
	        for(int i=1; i<length;i++) {
	        	
	        	if(min > recordList.get(i).getY()) {
	        		min = recordList.get(i).getY();
	        	}  	
	        }
	        
	        double regionYCenter = (min+max)/2;
	     /*   
	        System.out.println("::::::::::::::::::::::::::::");
	        System.out.println(regionXCenter);
	        System.out.println(regionYCenter);
	        System.out.println("::::::::::::::::::::::::::::");
	        */
	        for(Record record : recordList) {
	        	
					   minDistClusterArray[counter] = calculateDistance(record, regionXCenter, regionYCenter);
					   	x[counter] = record.getX();
					    y[counter] = record.getY();
					    nodeNumber[counter] = record.getNodeNumber();
					    counter++;
					    
				     }
	        
	        
	        int index=0;
	                
	        double smallest = minDistClusterArray[0];
	        
	        for(int i=0; i < counter; i++) {
	        	if(smallest > minDistClusterArray[i] ) {
	        		smallest = minDistClusterArray[i] ;
	        		index=i;
	        	}
	        } 

	        ClusterHeadSet clusterHead = new ClusterHeadSet(x[index], y[index], nodeNumber[index], key);
	        clusterHeadList.add(clusterHead);   
	        regionXCenter = regionXCenter+regionInterval; 
	      }
	}

	private double calculateDistance(Record record, double x, double y) {
		return Math.sqrt(Math.pow((x - record.getX()), 2) + Math.pow(y - record.getY(),2));
    }
	
	private void createClusterHeadAdjaceny() {
	
		int listLength;
		listLength = clusterHeadList.size();
		
		System.out.println("length"+listLength);
		
		double xArray[] = new double[listLength];
		double yArray[] = new double[listLength];
		double centroidArray[] = new double[listLength];;
		int i=0;
		
		for(ClusterHeadSet clusterHead : clusterHeadList) {
			 xArray[i] = clusterHead.getX();
			 yArray[i] = clusterHead.getY();
			 centroidArray[i] = clusterHead.getClusterHeadNumber();
			  i++;
		 }
		
		System.out.println(":::::::::::::::::::::::::::");
		for(int j=0;j<listLength;j++) {
			System.out.println( xArray[j]);
		}
		 printCombinationXArray(xArray, yArray, centroidArray, listLength, r);
		
	}
	
	public void printCombinationXArray(double xArray[], double yArray[], double centroidArray[], int listLength, int r) { 
		// A temporary array to store all combination one by one 
		double data[]=new double[r]; 
		double yData[] = new double[r];
		double centroidData[] = new double[r];
		
		// Print all combination using temprary array 'data[]' 
		combinationUtilXArray(xArray, yArray ,centroidArray, data, yData, centroidData, 0, listLength-1, 0, r); 
	} 

		
	public void combinationUtilXArray(double xArray[], double yArray[], double centroidArray[], double data[],  double yData[], double centroidData[] , int start,int end, int index, int r){ 
		if (index == r) 
		{ 
			double clusterHeadDistance;
			clusterHeadDistance = calculateEuclideanDistance(data, yData,r);
	
			adjacencyMatrix[ (int) centroidData[0]-1 ][(int) centroidData[1]-1 ] = clusterHeadDistance;
			adjacencyMatrix[(int) centroidData[1]-1 ][(int) centroidData[0]-1 ] = clusterHeadDistance;	
						
			return; 
		} 

		for (int i=start; i<=end && end-i+1 >= r-index;  i++) { 
			data[index] = xArray[i]; 
			yData[index] = yArray[i];
			centroidData[index] = centroidArray[i];
	        combinationUtilXArray(xArray, yArray, centroidArray, data, yData, centroidData, i+1, end, index+1, r); 
		} 
	} 
	
	public double calculateEuclideanDistance(double[] data, double[] yData, int r) {
		
		return Math.sqrt(Math.pow((data[0] -data[1]), 2) + Math.pow(yData[0] - yData[1],2));
	}
	
	private void printClusterHeadAdjaceny() {
		
		System.out.println(" Weighted Matrix  " );
		System.out.println(clusterNumber);
		
		for (int i=0; i<clusterNumber; i++) {
			for(int j=0; j<clusterNumber; j++) {
				System.out.print(adjacencyMatrix[i][j]+"  " );
			}
			
			System.out.println("  " );
			
		}
	}
	
	
	public void shortestPathAlgo() {
		
		double source, destination;
		double sourceToClusterHead = 0, destinationToclusterHead = 0;
		double pathDijkstra[] = new double[TOTAL_DATA-3];		
		
		Record sourceRecord = null, destinationRecord = null;
		ClusterHeadSet clusterHead = null;
		
		
		
		System.out.println("Enter the source node");
		//source = scan.nextDouble();
		source = 5;
		
		System.out.println("Enter the destination node");
		//destination = scan.nextDouble();
		destination =  9;
		
		for (Map.Entry<Integer, List<Record>> entry : clusterRecords.entrySet()) {
			
		     List<Record> recordList = entry.getValue();
		     
		     for(Record record : recordList) {
		    	 
		    	 
		    	 if(source == record.getNodeNumber()) {
		    		 sourceRecord = record;
		    		 System.out.println("sourceRecord"+sourceRecord);
		    	 }
		    	 if(destination == record.getNodeNumber()){
		    		
		    		 destinationRecord = record;
		    		 System.out.println("destinationRecord"+destinationRecord);
		    	 }
		    	 
		     }     
		     
		}
		
		int sCNumber, dCNumber;
		
		sCNumber = (int) sourceRecord.getClusterNumber();
		dCNumber = (int) destinationRecord.getClusterNumber();
		
		//System.out.println("sCNumber"+sCNumber+"dCNumber"+dCNumber);
		
		for(ClusterHeadSet clusterHeadSet : clusterHeadList) {
			
			if(sCNumber == clusterHeadSet.getClusterHeadNumber()){
				clusterHead = clusterHeadSet;
				sourceToClusterHead = calculateEuclideanDistanceCluster(clusterHead,sourceRecord);
				
			}
			else if(dCNumber == clusterHeadSet.getClusterHeadNumber()){
				destinationToclusterHead = calculateEuclideanDistanceCluster(clusterHeadSet,destinationRecord);
				
			}
		}
		
		//System.out.println(clusterHead);
		
		pathDijkstra[pathCounter] = sourceRecord.getNodeNumber();
		pathCounter++;
		
		//Source and destination within the same cluster
		if(sCNumber == dCNumber)
		{
			for(ClusterHeadSet clusterHeadSet : clusterHeadList) {
				if(clusterHeadSet.getClusterHeadNumber() == sourceRecord.getClusterNumber()) {
					pathDijkstra[pathCounter] = clusterHeadSet.getNodeNumber();
					pathCounter++;
					
					if(clusterHeadSet.getNodeNumber()!=destinationRecord.getNodeNumber()) {
						pathDijkstra[pathCounter] = destinationRecord.getNodeNumber();
						pathCounter++;
					}
					
				}
				
			}
			
			//pathDijkstra[pathCounter];
			
			totalDistance = sourceToClusterHead+destinationToclusterHead;
			outsideClusterDistance = totalDistance;
			insideClusterDistance = 0;
			
			System.out.println();
			System.out.println("Shortest Path from Source to Destination");
			for(int i=0;i<pathCounter;i++) {
				System.out.print((int) pathDijkstra[i]+" ");
			}
			System.out.println();
			
			
			System.out.println("Shorstest Path");
			System.out.println(totalDistance);
			
		}
		else {
			
	
		dijkstra(adjacencyMatrix, sCNumber-1, dCNumber-1, pathDijkstra);
		pathDijkstra[pathCounter] = destinationRecord.getNodeNumber();
		pathCounter++;
		
		totalDistance = sourceToClusterHead+destinationToclusterHead+shortestPathInsideCluster;
		insideClusterDistance = sourceToClusterHead+destinationToclusterHead;
		outsideClusterDistance = shortestPathInsideCluster;
		
		System.out.println();
		System.out.println("Shortest Path from Source to Destination");
		for(int i=0;i<pathCounter;i++) {
			System.out.print((int) pathDijkstra[i]+" ");
		}
		System.out.println();
		
		
		System.out.println("Shorstest Path");
		System.out.println(totalDistance);
		}
		
	}
	
	public double calculateEuclideanDistanceCluster(ClusterHeadSet clusterHead, Record record) {
		
		return Math.sqrt(Math.pow((record.getX() - clusterHead.getX()), 2) + Math.pow(record.getY() - clusterHead.getY(),2));
	}
	
	
	public void dijkstra(double[][] adjacencyMatrix, int startVertex,int dest, double[] pathDijkstra) { 
		int nVertices = adjacencyMatrix[0].length; 
		
		// shortestDistances[i] will hold the 
		// shortest distance from src to i 
		double[] shortestDistances = new double[nVertices]; 
		
		// added[i] will true if vertex i is 
		// included / in shortest path tree 
		// or shortest distance from src to  
		// i is finalized 
		boolean[] added = new boolean[nVertices]; 
		
		// Initialize all distances as  
		// INFINITE and added[] as false 
		for (int vertexIndex = 0; vertexIndex < nVertices; vertexIndex++) { 
		shortestDistances[vertexIndex] = Double.MAX_VALUE; 
		added[vertexIndex] = false; 
		} 
		
		// Distance of source vertex from itself is always 0 
		shortestDistances[startVertex] = 0; 
		
		// Parent array to store shortest path tree 
		double[] parents = new double[nVertices]; 
		
		// The starting vertex does not have a parent 
		parents[startVertex] = NO_PARENT; 
		
		// Find shortest path for all vertices 
		for (int i = 1; i < nVertices; i++) { 
		
		// Pick the minimum distance vertex from the set of vertices not yet processed.
		// nearestVertex is always equal to startNode in first iteration. 
		int nearestVertex = -1; 
		double shortestDistance = Double.MAX_VALUE; 
		for (int vertexIndex = 0; vertexIndex < nVertices; vertexIndex++) { 
		
			if (!added[vertexIndex] && shortestDistances[vertexIndex] <  shortestDistance)  
			{ 
				nearestVertex = vertexIndex; 
				shortestDistance = shortestDistances[vertexIndex]; 
			} 
		} 
		
		// Mark the picked vertex as processed 
		added[nearestVertex] = true; 
		
		// Update dist value of the adjacent vertices of the picked vertex. 
		for (int vertexIndex = 0; vertexIndex < nVertices; vertexIndex++)  { 
			double edgeDistance = adjacencyMatrix[nearestVertex][vertexIndex]; 
		
		if (edgeDistance > 0 && ((shortestDistance + edgeDistance) <  shortestDistances[vertexIndex]))  { 
			
			parents[vertexIndex] = nearestVertex; 
			shortestDistances[vertexIndex] = shortestDistance + edgeDistance; 
				} 
			} 
		} 
		
		printSolution(startVertex, shortestDistances, parents,dest,pathDijkstra); 
	} 
	
	// A utility function to print the constructed distances array and shortest paths 
	public void printSolution(int startVertex, double[] distances, double[] parents,int dest, double[] pathDijkstra) 
	{ 
		int nVertices = distances.length; 
		/*
		System.out.print("Vertex\t Distance\tPath"); 
		
		System.out.print("\n" + startVertex + " -> "); 
		System.out.print(dest + " \t\t "); 
		System.out.print(distances[dest] + "\t\t"); 
		*/
		shortestPathInsideCluster = distances[dest];
		printPath(dest, parents, pathDijkstra);
	} 
	
	// Function to print shortest path from source to currentVertex using parents array 
	public void printPath(double currentVertex, double[] parents, double[] pathDijkstra) { 
	
		// Base case : Source node has been processed 
		if (currentVertex == NO_PARENT) 
		{ 
		return; 
		} 
		printPath(parents[(int)currentVertex], parents, pathDijkstra); 
		
		System.out.print(currentVertex + " "); 
		
		for(ClusterHeadSet clusterHeadSet : clusterHeadList) {
			if(currentVertex+1 == clusterHeadSet.getClusterHeadNumber()) {
				pathDijkstra[pathCounter] = clusterHeadSet.getNodeNumber();
				pathCounter++;
			}
		}
		
		
	} 
	
	
	void printResult() {
		
		
		System.out.println();
		System.out.println();
		System.out.println("#################################################");
		
		double d=100;
		for(int l=0;l<5;l++) {
		System.out.println("Size of packet in bits");
		double packet = d;
		
		double nBW = Math.pow(10, 6);
		double BW = 5*nBW;
	
		double speed =  Math.pow(10, 6);
		double opticalSpeed = 2.1*Math.pow(10, 8);
		double transTime=0;
		
		int trackLinks;
		trackLinks = pathCounter-2;
		
		if(trackLinks==1) {
			transTime = (packet/nBW);
		}
		else {
	
			for(int j=0;j<2;j++) {
				transTime = transTime+(packet/nBW);
			}
			
			for(int i=0;i<trackLinks-2;i++) {
				transTime = transTime+(packet/BW);
			}
		}
		double propTime;
		propTime = (insideClusterDistance/opticalSpeed)+(outsideClusterDistance/speed);
		double efficiency;
		
		double queueDelay = 2*propTime;
		efficiency = transTime/(transTime+propTime+queueDelay);
		double m = efficiency*100;
		
		double throughput = efficiency*BW;
		double x = (transTime+propTime);	
		double y = throughput/Math.pow(10, 6);
		

		double a=(insideClusterDistance/opticalSpeed);
		double b=(outsideClusterDistance/speed);
		
		double g = 10*Math.pow(10, 6);
		double v = 20*Math.pow(10, 6);
		double c = 30*Math.pow(10, 6);
		double f = 40*Math.pow(10, 6);
		double e = 50*Math.pow(10, 6);
		
		
		double RTTa = 2*(totalDistance/g);
		double RTTb = 2*(totalDistance/v);
		double RTTc = 2*(totalDistance/c);
		double RTTf = 2*(totalDistance/f);
		double RTTe = 2*(totalDistance/e);
		
		
		System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11");
		System.out.println("RTTa="+ RTTa);
		System.out.println("RTTb="+RTTb);
		System.out.println("RTTc="+RTTc);
		System.out.println("RTTf="+RTTf);
		System.out.println("RTTe="+RTTe);
		System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		
		System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11");
		System.out.println("DATA="+ d);
		System.out.println("pathcounter"+pathCounter);
		System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		
		System.out.println("insideClusterDistance/opticalSpeed"+a);
		System.out.println("outsideClusterDistance/speed"+b);
		
		System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
		
		System.out.println("Transmission Time="+ transTime);
		System.out.println("Propogation Time="+ propTime);
		System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
		System.out.println("Efficiency="+efficiency);
		System.out.println("Throughput="+throughput);
		System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
		System.out.println("Efficiency in Percentage="+m);
		System.out.println("Throughput in Mbps="+y);
		System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
		System.out.println("RTT="+(propTime*2));
		System.out.println("Total Cycle Time of Algorithm for Data Transmission="+(transTime+propTime+(insideClusterDistance/opticalSpeed)));
		System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
		
		d=d+100;
		}
		
	}

	

}
