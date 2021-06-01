
public class Record {
	
	private double x;
	private double y;
	private double nodeNumber;
	private double clusterNumber;

	public Record(double x, double y, double nodeNumber) {
		super();
		this.x = x;
		this.y = y;
		this.nodeNumber = nodeNumber;
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

	public double getClusterNumber() {
		return clusterNumber;
	}
	public void setClusterNumber(double d) {
		this.clusterNumber = d;
	}
	
	@Override
	public String toString() {
		return "Record [x=" + x + ", y=" + y  + ", nodeNumber=" + nodeNumber + ",  clusterNumber="	+ clusterNumber + "]";
	}
	

}
