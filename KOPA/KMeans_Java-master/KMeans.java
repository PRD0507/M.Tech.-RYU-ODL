import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;


public class KMeans {

	private static final int NO_PARENT = -1; 
	static int clusterNumber = 3;					//update the cluster number based on the requirement
	List<Record> data = new ArrayList<Record>();
	List<Cluster> clusters = new ArrayList<Cluster>();
	List<ClusterHeadSet> clusterHeadList = new ArrayList<ClusterHeadSet>();
	double[][] adjacencyMatrix = new double[clusterNumber][clusterNumber] ;
	
		
	
	Map<Cluster, List<Record>> clusterRecords = new HashMap<Cluster, List<Record>>();
	int TOTAL_DATA;
	int r = 2;
	int pathCounter = 0;
	double shortestPathInsideCluster = 0;
	
	public static void main(String[] args) {
		long start = System.nanoTime();
		
		
		//adjacencyMatrix = new double[clusterNumber][clusterNumber];
		
		KMeans demo = new KMeans();
		demo.genereateRecord();
		demo.initiateClusterAndCentroid(clusterNumber);
		//demo.printRecordInformation();
		demo.printClusterInformation();
		demo.clusterHead();
		demo.printClusterHeadInformation();
		demo.createClusterHeadAdjaceny();
		demo.printClusterHeadAdjaceny();
		demo.shortestPathAlgo();
		
		long end = System.nanoTime();

		long execution = end - start;
		System.out.println("Execution time: " + execution + " nanoseconds");
		
	
	}


	private void genereateRecord() {
		//First three record mark the centroids
		//Add centroid record if you want to increase the clusters
		Record record = new Record(1.0, 1.0,0);
		data.add(record);
		record = new Record(2.5, 3.5,0);
		data.add(record);
		record = new Record(5.0, 7.0,0);
		data.add(record);
		//record = new Record(7.0,9.0,0);
		//data.add(record);
		
		//Record Data
		record = new Record(1.0, 1.0,1);
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
		record = new Record(8.0, 9.0,10);
		data.add(record);
		record = new Record(7.0, 10.0,11);
		data.add(record);
		
	}

	private void initiateClusterAndCentroid(int clusterNumber) {
		int counter = 1;
		Iterator<Record> iterator = data.iterator();
		//Iterator<Centroid> citerator = cData.iterator();
		TOTAL_DATA = data.size();
		
		
		Record record = null;
			
		
		while(iterator.hasNext()) {
			
			record = iterator.next();
			if(counter <= clusterNumber) {
				record.setClusterNumber(counter);
				initializeCluster(counter, record);
				counter++;
			
			}else {
				System.out.println(record);
				System.out.println("** Cluster Information **");
				for(Cluster cluster : clusters) {
					System.out.println(cluster);
				}
				System.out.println("*********************");
                double minDistance = Integer.MAX_VALUE;
                Cluster whichCluster = null;
                
                for(Cluster cluster : clusters) {
                	double distance = cluster.calculateDistance(record);
                	//System.out.println(distance);
                	if(minDistance > distance) {
                		minDistance = distance;
                		whichCluster = cluster;
                	}
                }
                
                record.setClusterNumber(whichCluster.getClusterNumber());
				whichCluster.updateCentroid(record);
				clusterRecords.get(whichCluster).add(record);

			}
			
			System.out.println("** Cluster Information **");
			for(Cluster cluster : clusters) {
				System.out.println(cluster);
			}
			System.out.println("*********************");
			
		}
	}

	private void initializeCluster(int clusterNumber, Record record) {
		
		Cluster cluster = new Cluster(clusterNumber,record.getX(),record.getY());
		clusters.add(cluster);
		List<Record> clusterRecord = new ArrayList<Record>();
		clusterRecord.add(record);
		clusterRecords.put(cluster, clusterRecord);

	}

	private void printRecordInformation() {
		   System.out.println("****** Each Record INFORMATIN *********");
		   for(Record record : data) {
			   System.out.println(record);
		   }
	   }

	private void printClusterInformation() {
	   System.out.println("****** FINAL CLUSTER INFORMATIN *********");
	   for (Map.Entry<Cluster, List<Record>> entry : clusterRecords.entrySet())  {
        System.out.println("Key = " + entry.getKey() + 
                         ", Value = " + entry.getValue()); 
	   }
	}
	
	private void printClusterHeadInformation() {
		   System.out.println("****** Each Cluster Head INFORMATION *********");
		   for(ClusterHeadSet clusterHead : clusterHeadList) {
			 //  System.out.println(clusterHead.getX());
			//   System.out.println(clusterHead.getY());
			//   System.out.println(clusterHead.getClusterHeadNumber());
			   System.out.println(clusterHead);
		   }
	   }
	

	private void clusterHead() {
	
		double xcentroid, ycentroid, clusterNumberArray[] = null;
		double clusterNumber;
		int flag=0;
		
		
		for (Map.Entry<Cluster, List<Record>> entry : clusterRecords.entrySet())  {
			
			
			int counter=0;
			double minDistClusterArray[] = new double[TOTAL_DATA];
			double x[] = new double[TOTAL_DATA];
			double y[] = new double[TOTAL_DATA];
			double nodeNumber[] = new double[TOTAL_DATA];
			
	        Cluster key = entry.getKey();
	        List<Record> recordList = entry.getValue();
	        
	      
	       
	        
	        flag=0;
	        for(Record record : recordList) {
	        	//Skipping the first centroid record
				  // System.out.println(record);
				   if(flag == 0) {
					   flag=1;
				
				   }
				   else {
					   minDistClusterArray[counter] = key.calculateDistance(record);
					  //  System.out.println("min"+ key.calculateDistance(record));
					   	x[counter] = record.getX();
					    y[counter] = record.getY();
					    nodeNumber[counter] = record.getNodeNumber();
					    counter++;
				   }
				  
			   }
	        
	        
	        int index=0;
	                
	        double smallest = minDistClusterArray[0];
	        
	        /*double element : minDistClusterArray*/
	        for(int i=0; i < counter; i++) {
	        	if(smallest > minDistClusterArray[i] ) {
	        		smallest = minDistClusterArray[i] ;
	        		index=i;
	        		
	        		//System.out.println("mind is"+ smallest);
	        	}
	        } 

	     
	        
	        ClusterHeadSet clusterHead = new ClusterHeadSet(x[index], y[index], nodeNumber[index], (int) key.getClusterNumber());
	        clusterHeadList.add(clusterHead); 
	      }
	}


	
	private void createClusterHeadAdjaceny() {
	
		int listLength;
		
		listLength = clusterHeadList.size();
		
		//System.out.println(clusterHeadList);		
		//printCombination(clusterHeadList, listLength, r); 
		
		
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
		 
		// printCombination(clusterHeadList, listLength, r);
		
	}
	
	
	public void printCombinationXArray(double xArray[], double yArray[], double centroidArray[], int listLength, int r) 
	{ 
	// A temporary array to store all combination one by one 
	double data[]=new double[r]; 
	double yData[] = new double[r];
	double centroidData[] = new double[r];
	
	// Print all combination using temprary array 'data[]' 
	combinationUtilXArray(xArray, yArray ,centroidArray, data, yData, centroidData, 0, listLength-1, 0, r); 
	} 
	
	public void printCombination(List<ClusterHeadSet> clusterHeadList, int listLength, int r) 
	{ 
	// A temporary array to store all combination one by one 
	List<ClusterHeadSet> clusterHeadListData =new ArrayList<ClusterHeadSet>();
	
	// Print all combination using temprary array 'data[]' 
	combinationUtil(clusterHeadList, clusterHeadListData, 0, listLength-1, 0, r); 
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
		for (int i=0; i<clusterNumber; i++) {
			for(int j=0; j<clusterNumber; j++) {
				System.out.print(adjacencyMatrix[i][j]+"  " );
			}
			
			System.out.println("  " );
			
		}
	}

	
	public void combinationUtil(List<ClusterHeadSet> clusterHeadList, List<ClusterHeadSet> clusterHeadListRecord, int start,int end, int index, int r){ 
		if (index == r) 
		{ 
			for (int j=0; j<r; j++) 
			System.out.print(clusterHeadListRecord.get(j)+" sada"); 
			System.out.println(""); 
			return; 
		} 

	for (int i=start; i<=end && end-i+1 >= r-index;  i++) { 
		
		/*
		double x, y, clusterNumber;
		x = clusterHeadList.get(index).getX();
		y = clusterHeadList.get(index).getY();
		clusterNumber = clusterHeadList.get(index).getClusterHeadNumber();
		*/
		
		clusterHeadListRecord.add(index, clusterHeadList.get(index));
		
	//	System.out.println("New infor");
	//	System.out.println(clusterHeadListRecord);
		combinationUtil(clusterHeadList, clusterHeadListRecord, i+1, end, index+1, r); 
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
		
		  
		for (Map.Entry<Cluster, List<Record>> entry : clusterRecords.entrySet()) {
			
			 Cluster key = entry.getKey();
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
			if(dCNumber == clusterHeadSet.getClusterHeadNumber()){
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
		
		System.out.print("Vertex\t Distance\tPath"); 
		
		System.out.print("\n" + startVertex + " -> "); 
		System.out.print(dest + " \t\t "); 
		System.out.print(distances[dest] + "\t\t"); 
		
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
		
		System.out.print(currentVertex+1 + " "); 
		
		for(ClusterHeadSet clusterHeadSet : clusterHeadList) {
			if(currentVertex+1 == clusterHeadSet.getClusterHeadNumber()) {
				pathDijkstra[pathCounter] = clusterHeadSet.getNodeNumber();
				pathCounter++;
			}
		}
		
		
		
		
	} 
	
	
}
		



	

	
	

	
	



	
	




