
public class Cluster {

//	private double xCentroid;
//	private double yCentroid;

	private double clusterNumber;
	
	

	public Cluster(double clusterNumber) {
		super();
		this.clusterNumber = clusterNumber;
		
	}
	

	

	public double getClusterNumber() {
		return clusterNumber;
	}


	public void setClusterNumber(double clusterNumber) {
		this.clusterNumber = clusterNumber;
	}


	@Override
	public String toString() {
		return "Cluster [ clusterNumber=" + clusterNumber + "]";
	}
		

}
