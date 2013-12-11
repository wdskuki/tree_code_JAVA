package treeCode;

//import first.Node;

public class Node  implements Cloneable{

	/**
	 * @param var
	 */
	

	//public int flowPatternsNum;
	public int pi; //�ýڵ���tree�ĸ��ڵ�
	public double flow;//��tree�иýڵ㵽���ڵ������
	int ID;
	
	/* pi��flow��ӦԪ�أ�
	 *  0: star
	 *  1: singleMaxSpanningTree
	 *  2: singleMinTimeTreeHeuristicFromMaxST
	 *  3: singleMinTimeTreeHeuristicFromStar
	 *  4: singleMinTimeTreeHeuristicIncludeNodeOneByOne
	 *  
	 *  5: singleMinTimeTreeHeuristicFromMaxSTchgFlow
	 *  6: singleMinTimeTreeHeuristicFromStarchgFlow
	 *  7:singleMinTimeTreeHeuristicIncludeNodeOneByOne
	 *  8: Jun's singleMaxSpanningTree strategy
	 */
	

	
	public Node(int ID){
		this.ID = ID;
	}
	protected Object clone(){
		// TODO Auto-generated method stub
		Node o = null;
		try{
			o = (Node)super.clone();
		}catch(CloneNotSupportedException e){
			e.printStackTrace();
		}
		return o;
	}
	public void Print(){
		System.out.print(ID + "\t");
		System.out.print(pi + "\t" + (short)flow + "\t\t");
		System.out.print("\n");
	}
	public static void main(String[] args) {
		// test: print node information
		int nodeNum = 10;
	
		Node[] node = new Node[nodeNum];
		for(int i = 0; i < nodeNum; i++)
			node[i] = new Node(i);
		System.out.print("ID" + "\t" );
		System.out.print("pi" + "\t" + "flow" + "\t\t");
		System.out.print("\n");
		for(int i = 0; i < nodeNum; i++)
			node[i].Print();
	}
}
