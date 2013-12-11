package treeCode;

public class Link {

	/**
	 * @param args
	 */
	public int point1; 
	public int point2;
	public double bandwidth;//bottleneckTime对应的带宽
	public double bottleneckTime;
	public double bottleneckBW;//tree里的最小带宽
	
	public Link(){
		point1 = -1;
		point2 = -1;
		bandwidth = 0;
		bottleneckTime = 0;
		bottleneckBW = 0;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
