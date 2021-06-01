
public class Cluster {

	private double xCentroid;
	private double yCentroid;

	private double clusterNumber;
	
	
	
	public Cluster(double clusterNumber, double xCentroid, double yCentroid) {
		super();
		this.clusterNumber = clusterNumber;
		this.xCentroid = xCentroid;
		this.yCentroid = yCentroid;
	}
	

	public double getxCentroid() {
		return xCentroid;
	}


	public void setxCentroid(double d) {
		this.xCentroid = d;
	}


	public double getyCentroid() {
		return yCentroid;
	}


	public void setyCentroid(double d) {
		this.yCentroid = d;
	}


	public double getClusterNumber() {
		return clusterNumber;
	}


	public void setClusterNumber(double clusterNumber) {
		this.clusterNumber = clusterNumber;
	}


	@Override
	public String toString() {
		return "Cluster [xCentroid=" + xCentroid + ", yCentroid=" + yCentroid + ", clusterNumber=" + clusterNumber + "]";
	}
	
	// Euclidean distance calculation
	public double calculateDistance(Record record) {
		return Math.sqrt(Math.pow((getxCentroid() - record.getX()), 2) + Math.pow(getyCentroid() - record.getY(),2));
    }

	// K-Mean Algorithm
	public void updateCentroid(Record record) {
		setxCentroid((getxCentroid()+record.getX())/2);
		setyCentroid((getyCentroid()+record.getY())/2);
		
	}
	
	

}
