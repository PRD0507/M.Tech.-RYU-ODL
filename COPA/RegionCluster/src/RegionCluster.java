import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;


public class RegionCluster {

	private static final int NO_PARENT = -1; 
	static int clusterNumber = 4;			//Update the cluster number based on the requirement
	int TOTAL_DATA;
	int r = 2;
	double minX, minY;
	double maxX, maxY;
	double xLineDistance;
	double regionInterval;
	int pathCounter = 0;
	double shortestPathInsideCluster = 0;
	double bandwidth;
	
	List<Record> data = new ArrayList<Record>();
	List<Cluster> clusters = new ArrayList<Cluster>();
	List<ClusterHeadSet> clusterHeadList = new ArrayList<ClusterHeadSet>();
	double[][] adjacencyMatrix = new double[clusterNumber][clusterNumber] ;
	
	

	Map<Integer, List<Record>> clusterRecords = new HashMap<Integer, List<Record>>();
	
	public static void main(String args[]) {
		
		RegionCluster demo = new RegionCluster();
		demo.genereateRecord();
		demo.initiateClusterAndCentroid(clusterNumber);
		demo.printClusterInformation();
		
		demo.clusterHead();
		demo.printClusterHeadInformation();
		demo.createClusterHeadAdjaceny();
		demo.printClusterHeadAdjaceny();
		demo.shortestPathAlgo();
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
		Record record = new Record(0.0, 1.0,1);
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
		data.add(record);
		
	}
	
	private void initiateClusterAndCentroid(int clusterNumber) {
			
		Iterator<Record> iterator = data.iterator();
		//Iterator<Centroid> citerator = cData.iterator();
		TOTAL_DATA = data.size();
		int counter = 1;
		
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
		regionInterval = xLineDistance/clusterNumber;
		
		//System.out.println(xLineDistance);
		//System.out.println(regionInterval);
		
		int len = data.size();
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
			/*	
			System.out.println("////////////////////////");
			System.out.println("data");
			
			int len = data.length;
			for(int m=0;m<len;m++) {
				System.out.println(data[m]);
			}
			System.out.println("ydata");
			
			len=yData.length;
			for(int m=0;m<len;m++) {
				System.out.println(yData[m]);
			}
			System.out.println("distance"+clusterHeadDistance);*/
			
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
		double totalDistance;
		double pathDijkstra[] = new double[TOTAL_DATA-3];		
		
		Record sourceRecord = null, destinationRecord = null;
		ClusterHeadSet clusterHead = null;
		
		Scanner scan = new Scanner(System.in);
		
		System.out.println("Enter the source node");
		source = scan.nextDouble();
		
		System.out.println("Enter the destination node");
		destination = scan.nextDouble();
		
		/*
		System.out.println("Enter the link bandwidth");
		bandwidth = scan.nextDouble();*/
		  
		for (Map.Entry<Integer, List<Record>> entry : clusterRecords.entrySet()) {
			
			 Integer key = entry.getKey();
		     List<Record> recordList = entry.getValue();
		     
		     for(Record record : recordList) {
		    	 if(source == record.getNodeNumber()) {
		    		 sourceRecord = record;
		    	 }
		    	 else if(destination == record.getNodeNumber()){	
		    		 destinationRecord = record;
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
	

	

}
