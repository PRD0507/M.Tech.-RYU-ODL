
public class ClusterHeadSet {

	private double x;
	private double y;
	private double nodeNumber;
	private double ClusterHeadNumber;
	
	public ClusterHeadSet(double x, double y, double nodeNumber, double clusterNumber) {
		super();
		this.x = x;
		this.y = y;
		this.nodeNumber = nodeNumber;
		this.ClusterHeadNumber = clusterNumber;
	}
	
	public double getNodeNumber() {
		return nodeNumber;
	}

	public void setNodeNumber(double nodeNumber) {
		this.nodeNumber = nodeNumber;
	}

	public double getX() {
		return x;
	}
	public void setX(double x) {
		this.x = x;
	}
	public double getY() {
		return y;
	}
	public void setY(double y) {
		this.y = y;
	}
	public double getClusterHeadNumber() {
		return ClusterHeadNumber;
	}
	public void setClusterHeadNumber(double ClusterHeadNumber) {
		this.ClusterHeadNumber = ClusterHeadNumber;
	}

	@Override
	public String toString() {
		return "ClusterHeadSet [x=" + x + ", y=" + y + ", NodeNumber= " + nodeNumber+ ", ClusterHeadNumber=" + ClusterHeadNumber + "]";
	}
	
}
