package treeCode;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Random;
import java.util.Map.Entry;

import FRC2.FRC2;

import treeCode.Link;
import treeCode.CodingPattern;
import treeCode.Node;

public class Tree {

	/**
	 * @param args
	 */
	//parameters in Regenerating Codes
	public int n;
	public int d;
	public int k;
	public double beta;
	public double alpha;
	public double M;	
	
	
	public double [][]W; //整个网络的带宽分布，是一个n*n的矩阵
	public int[] involvedNodeSet; //参与构造这棵树的节点集合（参与修复的节点集合 + 新节点（共有d+1个节点））
	public double[][] w;	//参与修复的节点集合 + 新节点（共有d+1个节点）的一个带宽分布，是一个nodeNum*nodeNum的矩阵，是W的一个子集合
	public Node[] V; //参与修复的节点集合 + 新节点（共有d+1个节点）
	public int nodeNum; //d+1,即参与节点数目+新节点
	public CodingPattern codingPattern; //MBR or MSR
	public int r; //involvedNodeSet[r]表示new comer
	//public int flowPatternNum;  
	public static final int MAX_VALUE = 10000000;
	
	
	//一些基本参数的初始化
 	public Tree(int n, int d, int k, double M, CodingPattern codingPattern,  int r){
		this.n = n;
		this.d = d;
		this.k = k;
		this.M = M;
		this.codingPattern = codingPattern;
		//this.flowPatternNum = flowPatternNum;
		switch(codingPattern){
			case MSR: 	Msr();break;
			case MBR: Mbr();break;
			default: System.out.println("wrong Coding Pattern!");
		}
		
		this.r = r;
		nodeNum = d+1;
		W = new double[n][n];
		w = new double[nodeNum][nodeNum];
		involvedNodeSet = new int[nodeNum];
		V = new Node[nodeNum];
		for(int i = 0; i < nodeNum; i++)
			V[i] = new Node(i);
	}
	
 	//所有节点之间两两之间的带宽
 	public void set_W(double [][]W){
 		for(int i = 0; i < n; i++){
 			for(int j = 0; j < n; j++){
 				this.W[i][j] = W[i][j];
 			}
 		}
 	}
 	
 	//设置new comer r
 	public void set_r(int r){
 		this.r = r;
 	}
 	
 	//设置参与节点的编号
 	public void set_involvedNodeSet(int[] provider_node_id){
 		for(int i = 0; i < nodeNum; i++){
 			this.involvedNodeSet[i] = provider_node_id[i];
 		}
 	}
 	
 	//获得新节点和参与节点之间的两两带宽值
 	public void get_w_according_to_involvedNodeSet(){
 		for(int i = 0; i < nodeNum; i++){
 			for(int j = 0; j < nodeNum; j++){
 				int i1 = involvedNodeSet[i];
 				int j1 = involvedNodeSet[j];
 				w[i][j] = W[i1][j1];
 			}
 		}
 	}
 	
 	
 	/*test or debug*/
 	//从0~n-1个中选出d个providerNodes和1个failed node
 	 	public void random_set_providerNode_and_newNode(){
 	 		int[] seed = new int[n];
 	 		for(int i = 0; i < n; i++) seed[i] = i;
 	 		int len = seed.length;
 	 		
 	 		int[] result = new int[len];
 	 		Random random = new Random();
 	 		for(int i = 0; i < len; i++){
 	 			int r = random.nextInt(len - i);
 	 			result[i] = seed[r];
 	 			seed[r] = seed[len - 1 - i];
 	 		}
 	 		//System.out.print("result:" + Arrays.toString(result));
 	 		for(int i = 0; i < nodeNum; i++)  involvedNodeSet[i] = result[i];
 	 		set_r(0); //involveNodeSet集合中的第一个节点作为failed节点
 	 //		System.out.println("involvedNodeSet:" + Arrays.toString(involvedNodeSet));
 	 	}
 		public void randomDirectedLinkbandwidth_W(double l, double h){
 			double low = l;
 			double high = h;
 			Random r = new Random();
 			for(int i = 0; i < n; i++){
 				for(int j = 0; j < n; j++){
 					W[i][j] = r.nextDouble()*(high - low) + low;
 				}
 				W[i][i] = -1;
 			}
 		}
 		public void randomIndirectedLinkBandwidth_W(double l, double h){
 			double low = l;
 			double high = h;
 			Random r = new Random();
 			for(int i = 0; i < n; i++){
 				for(int j = i+1; j < n; j++){
 					W[i][j] = W[j][i] = r.nextDouble()*(high - low) + low;
 				}
 				W[i][i] = -1;
 			}
 		}
 		public void randomDirectedLinkbandwidth_W(double[] L, double[] H, int num){
 			Random r1 = new Random();
 			Random r2 = new Random();
 			for(int i = 0; i < n; i++){
 				for(int j = 0; j < n; j++){
 					int r11 = r1.nextInt(num);
 					double r22 = r2.nextDouble()*(H[r11] - L[r11]) + L[r11];
 					W[i][j] = r22;
 				}
 			}
 		}
 		public void randomIndirectedLinkbandwidth_W(double[] L, double[] H, int num){
 			Random r1 = new Random();
 			Random r2 = new Random();
 			for(int i = 0; i < n; i++){
 				for(int j = i+1; j < n; j++){
 					int r11 = r1.nextInt(num);
 					double r22 = r2.nextDouble()*(H[r11] - L[r11]) + L[r11];
 					W[i][j] =W[j][i] =  r22;
 				}
 			}
 		}
 		
 		public void testIndirectedLinkBandwidth(){
 			nodeNum = 5;
 			w[0][1] = w[1][0] = 40;
 			w[0][2] = w[2][0] = 20;
 			w[0][3] = w[3][0] = 3;	//
 			w[0][4] = w[4][0] = 30;
 			w[1][2] = w[2][1] = 15;
 			w[1][3] = w[3][1] = 45;
 			w[1][4] = w[4][1] = 50;
 			w[2][3] = w[3][2] = 10;
 			w[2][4] = w[4][2] = 40;
 			w[3][4] = w[4][3] = 35;
 		}
 		public void testIndirectedLinkBandwidth2(){
 			nodeNum = 5;
 			w[0][1] = w[1][0] = 40;
 			w[0][2] = w[2][0] = 20;
 			w[0][3] = w[3][0] = 25;	//
 			w[0][4] = w[4][0] = 30;
 			w[1][2] = w[2][1] = 15;
 			w[1][3] = w[3][1] = 45;
 			w[1][4] = w[4][1] = 50;
 			w[2][3] = w[3][2] = 10;
 			w[2][4] = w[4][2] = 40;
 			w[3][4] = w[4][3] = 35;
 		}
 		public void testIndirectedLinkBandwidth3(){
 			nodeNum = 5;
 			w[0][1] = w[1][0] = 80;
 			w[0][2] = w[2][0] = 19;
 			w[0][3] = w[3][0] = 25;	//
 			w[0][4] = w[4][0] = 40;
 			w[1][2] = w[2][1] = 21;
 			w[1][3] = w[3][1] = 45;
 			w[1][4] = w[4][1] = 39;
 			w[2][3] = w[3][2] = 10;
 			w[2][4] = w[4][2] = 40;
 			w[3][4] = w[4][3] = 30;
 		}	
 		public void testIndirectedLinkBandwidth4(){
 			nodeNum = 5;
 			w[0][1] = w[1][0] = 41;
 			w[0][2] = w[2][0] = 21;
 			w[0][3] = w[3][0] = 25;	//
 			w[0][4] = w[4][0] = 40;
 			w[1][2] = w[2][1] = 21;
 			w[1][3] = w[3][1] = 41;
 			w[1][4] = w[4][1] = 39;
 			w[2][3] = w[3][2] = 10;
 			w[2][4] = w[4][2] = 40;
 			w[3][4] = w[4][3] = 50;
 		}	
 		public void testIndirectedLinkBandwidth5(){
 			nodeNum = 5;
 			w[0][1] = w[1][0] = 6;
 			w[0][2] = w[2][0] = 2;
 			w[0][3] = w[3][0] = 2.5;	//
 			w[0][4] = w[4][0] = 3;
 			w[1][2] = w[2][1] = 15;
 			w[1][3] = w[3][1] = 45;
 			w[1][4] = w[4][1] = 50;
 			w[2][3] = w[3][2] = 10;
 			w[2][4] = w[4][2] = 40;
 			w[3][4] = w[4][3] = 35;
 		}
 		public void testIndirectedLinkBandwidth6(){
 			nodeNum = 5;
 			w[0][1] = w[1][0] = 10;
 			w[0][2] = w[2][0] = 2;
 			w[0][3] = w[3][0] = 2.5;	//
 			w[0][4] = w[4][0] = 3;
 			w[1][2] = w[2][1] = 15;
 			w[1][3] = w[3][1] = 45;
 			w[1][4] = w[4][1] = 50;
 			w[2][3] = w[3][2] = 10;
 			w[2][4] = w[4][2] = 40;
 			w[3][4] = w[4][3] = 35;
 		}
 		
 		
 	//star structure	
 	public Link starMinTime(){
		for(int i = 0; i < nodeNum; i++){
			V[i].pi = r;
			V[i].flow = beta;
		}
		V[r].pi = -1;
		V[r].flow = 0;
		Link bottleneckLink = bottleneckBandwidthLink();
		V[r].flow = 0;
		bottleneckLink.bottleneckBW = getTreeBottleneckBW();
		return bottleneckLink;
	}
	
 	//star structure with change flow in links
	public Link starMinTimeChgFlow(){
		starMinTime();
		clearLink();
		CompEachNodeFlow();
		V[r].pi = -1;
		V[r].flow = 0;
		Link bottleneckLink = bottleneckBandwidthLink();
		V[r].flow = 0;
		bottleneckLink.bottleneckBW = getTreeBottleneckBW();
		return bottleneckLink;
	}
 	
	//Tree structure
	public Link singleMinTimeTreeHeuristicIncludeNodeOneByOne(){
		Link bottleneckLink = new Link();
		LinkedList<Integer> A = new LinkedList<Integer>(); //the  candidate nodes of the tree(still not in the tree, but maybe in the tree)
		LinkedList<Integer> B = new LinkedList<Integer>(); //the nodes set in the tree
		V[r].flow = beta;
		V[r].pi = -1;
		B.add(Integer.valueOf(r));//add  newcomer in B
		for(int i = 0; i < nodeNum; i++){
			if(i != r)
				A.add(Integer.valueOf(i)); //add others nodes in A
		}
		while(A.size() > 0){
			int m_u = A.getFirst().intValue();
			int m_v = B.getFirst().intValue();
			double minTime = CompRouter4(m_u,m_v); //computer the min time from m_v to m_u
			
			//find the best nodes in B which is the father node(in B) of m_u
			for(int i = 0; i < A.size(); i++){
				int u = A.get(i).intValue();
				for(int j = 0; j < B.size(); j++){
					int v = B.get(j).intValue();
					double time = CompRouter4(u,v);
					if(minTime > time){
						minTime = time;
						m_u = u;
						m_v = v;
					}
				}
			}
			
			A.remove(Integer.valueOf(m_u));
			B.add(Integer.valueOf(m_u));
			V[m_u].pi = m_v;
			//System.out.println("n = " + n);
			V[m_u].flow = beta;
			
			//change all the links flow iteratively from m_u to root r, and find the bottleneckLink incidentally. 
			while(m_u != r){
				int x = V[m_u].pi;
				if(V[m_u].flow + V[x].flow <= alpha)
					V[x].flow = V[m_u].flow + V[x].flow;
				else
					V[x].flow = alpha;
				if(bottleneckLink.bottleneckTime < minTime && minTime == V[m_u].flow/w[m_u][x]){
					bottleneckLink.bottleneckTime = minTime;
					bottleneckLink.point1 = m_u;
					bottleneckLink.point2 = x;
					bottleneckLink.bandwidth = w[m_u][x];
				}
				m_u = x;
			}
		}
		V[r].pi = -1;
		V[r].flow = 0;
		bottleneckLink.bottleneckBW = getTreeBottleneckBW();
		return bottleneckLink;
	}
	
	//Tree structure with change flow in links
	public Link singleMinTimeTreeHeuristicIncludeNodeOneByOneChgFlow(){
		//CompEachNodeFlow在原有的singleMinTimeTreeHeuristicIncludeNodeOneByOne的基础上进行flow的调整
		CompEachNodeFlow();//flowPattern对应不同的启发式算法，详细可以参考Node的定义，这里为7，对应的就是singleMinTimeTreeHeuristicIncludeNodeOneByOneChgFlow
		V[r].pi = -1;
		V[r].flow = 0;
		Link bottleneckLink = bottleneckBandwidthLink();
		V[r].flow = 0;
		bottleneckLink.bottleneckBW = getTreeBottleneckBW();
		return bottleneckLink;
	}
	
	//Tree structure in Jun's Paper in IWQoS'09
	public Link JunSingleMaxSpanningTree(){
		generalMaxSpanningTree_Prim();
		for(int i = 0; i < nodeNum; i++)
			V[i].flow = beta;
		V[r].pi = -1;
		V[r].flow = 0;
		JunAdjustRootLink(d-k+1);
		//Link bottleneckLink = bottleneckBandwidthLink(flowPattern);
		Link bottleneckLink = new Link();
		bottleneckLink.bottleneckTime = 0;
		for(int i = 0; i < nodeNum; i++){
			if(i != r){
				int j = V[i].pi;
				double time = V[i].flow/w[i][j];
				if(bottleneckLink.bottleneckTime < time){
					bottleneckLink.bottleneckTime = time;
					bottleneckLink.point1 = i;
					bottleneckLink.point2 = j;
					bottleneckLink.bandwidth = w[i][j];
				}
			}
		}
		//V[r].flow[flowPattern] = 0;
		bottleneckLink.bottleneckBW = getTreeBottleneckBW();
		return bottleneckLink;
	}
	
	
	public double getTreeBottleneckBW(){
		double bottleneckBW = MAX_VALUE;
		for(int i = 0; i < nodeNum; i++){
			if(i != r){
				int j = V[i].pi;
				if(bottleneckBW > w[i][j])
					bottleneckBW = w[i][j];
			}
		}
		return bottleneckBW;
	}
	
	
	/*print function*/
	public void Print_bottleneck(Link bottleneckLink){
		System.out.println("------");
		System.out.println("bottleLinkPoint1 = " + bottleneckLink.point1);
		System.out.println("bottleLinkPoint2 = " + bottleneckLink.point2);
		System.out.println("bottleneckTime = " + bottleneckLink.bottleneckTime);
		System.out.println("bandwidth = " + bottleneckLink.bandwidth);
	}
	public void Print_Node(){
		System.out.println("------");
		System.out.print("ID" + "\t" );
		System.out.print("pi" + "\t" + "flow" + "\t\t");
		System.out.print("\n");
		for(int i = 0; i < nodeNum; i++)
			V[i].Print();
	}
 	
	private void Mbr(){
		alpha = 2*M*d/(k*(2*d - k + 1));
		beta = 2*M/(k*(2*d - k + 1));
	}
	private void Msr(){
		alpha = M/k;
		beta = M/(k*(d-k+1));
	}
	
	private Link bottleneckBandwidthLink(){
		Link bottleneckLink  = new Link();
		int[][] linkGraph = new int[nodeNum][nodeNum];
		for(int i = 0; i < nodeNum; i++){
			for(int j = 0; j < nodeNum; j++){
				linkGraph[i][j] = 0;
				if(i == V[j].pi){
					linkGraph[i][j] = 1;
				}
			}
		}
		
		outputFlow(r, linkGraph, bottleneckLink);
		
		return bottleneckLink;
//		double[] outputFlow = new double[nodeNum];
//		for(int i = 0; i < nodeNum; i++) 
//			outputFlow[i] = V[i].flow[flowPattern];
//		
//		for(int i = 0; i < nodeNum; i++){
//			int j = V[i].pi[flowPattern];
//			
//		}
		
	}
	private double outputFlow(int i, int[][] linkGraph, Link bottleneckLink){
		double totalFlow = 0.0;
		for(int j = 0; j < nodeNum; j++){
			if(linkGraph[i][j] == 1){
				totalFlow += outputFlow(j, linkGraph, bottleneckLink);
			}
		}
		totalFlow += V[i].flow;
		if(totalFlow > alpha)
			totalFlow = alpha;
		if(i != r){
			double time = totalFlow/w[i][V[i].pi];
			if(bottleneckLink.bottleneckTime < time){
				bottleneckLink.bottleneckTime = time;
				bottleneckLink.bandwidth = w[i][V[i].pi];
				bottleneckLink.point1 = i;
				bottleneckLink.point2 = V[i].pi;
			}
		}
		V[i].flow = totalFlow;
		return totalFlow;
	}
	

	private void CompEachNodeFlow(){		
		double[] bottleneckBandwidth = new double[nodeNum];
		double[] bottleneckBandwidth2 = new double[nodeNum];
		int[] flowNum = new int[nodeNum];
		for(int i = 0; i < nodeNum; i++){
			bottleneckBandwidth[i] = 0;
			bottleneckBandwidth2[i] = 0;
			flowNum[i] = 0;
		}
		for(int i = 0; i < nodeNum; i++){
			int x = i;
			//flowNum[x] = 1;
			while(x != r){
				flowNum[x] += 1;
				x = V[x].pi;
			}
		}
		for(int i = 0; i < nodeNum; i++){
			if(i != r){
				bottleneckBandwidth[i] = w[i][V[i].pi]/flowNum[i];
			}
		}

		
		for(int i = 0; i < nodeNum; i++){
			if(i != r){
				int x = i;
				bottleneckBandwidth2[i] = bottleneckBandwidth[x];
				while(x != r){
					if(bottleneckBandwidth2[i] > bottleneckBandwidth[x])
						bottleneckBandwidth2[i] = bottleneckBandwidth[x];
					x = V[x].pi;
				}
			}
		}
		bottleneckBandwidth2[r] = MAX_VALUE;
		//return bottleneckBandwidth2;
		if(codingPattern == CodingPattern.MBR)
			assignFlow2(bottleneckBandwidth2);
		if(codingPattern == CodingPattern.MSR)
			assignFlow(bottleneckBandwidth2);
	}
	
	
	private void assignFlow(double[] bottleneckBandwidth2){
		HashMap<Integer, Double> map = new HashMap<Integer, Double>();
		for(int i = 0; i < nodeNum; i++)
			map.put(Integer.valueOf(i), Double.valueOf(bottleneckBandwidth2[i]));

		ArrayList<Entry<Integer, Double> > al = new ArrayList<Entry<Integer, Double> > (map.entrySet());
		Collections.sort(al, new Comparator<Map.Entry<Integer, Double> >(){
			public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2){
				if (o2.getValue().doubleValue() > o1.getValue().doubleValue())
					return -1;
				else if(o2.getValue().doubleValue() == o1.getValue().doubleValue())
					return 0;
				else
					return 1;
			}
		});
		//System.out.println(al);
		
		double bandwidthSum = 0;
		for(int i = 0; i < d-k+1; i++)
			bandwidthSum += al.get(i).getValue();
		for(int i = 0; i < d-k+1; i++)
			V[al.get(i).getKey()].flow = alpha*(al.get(i).getValue()/bandwidthSum);
		for(int i = d-k+1; i < d; i++)
			V[al.get(i).getKey()].flow = V[al.get(i-1).getKey()].flow;
	}
	

	private void assignFlow2(double[] bottleneckBandwidth2){
		HashMap<Integer, Double> map = new HashMap<Integer, Double>();
		for(int i = 0; i < nodeNum; i++)
			map.put(Integer.valueOf(i), Double.valueOf(bottleneckBandwidth2[i]));

		ArrayList<Entry<Integer, Double> > al = new ArrayList<Entry<Integer, Double> > (map.entrySet());
		Collections.sort(al, new Comparator<Map.Entry<Integer, Double> >(){
			public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2){
				if (o2.getValue().doubleValue() > o1.getValue().doubleValue())
					return -1;
				else if(o2.getValue().doubleValue() == o1.getValue().doubleValue())
					return 0;
				else
					return 1;
			}
		});
		//System.out.println(al);
		double[] BW = new double[d];
		for(int i = 0; i < d; i++) BW[i] = al.get(i).getValue().doubleValue();
		
		Object[] result = null;
		try{
			FRC2 frc2 = new FRC2();
			result = frc2.NonUniformFlow(2, k, d, alpha, beta, BW);

			double lamda = Double.parseDouble(result[0].toString());
//			System.out.println(result[0]);
//			System.out.println(result[1]);

			String[] str2 = result[1].toString().split("\n");
//			System.out.println(str2[2]);
			if(codingPattern == CodingPattern.MBR){
				for(int i = 0; i < d; i++)
					V[al.get(i).getKey()].flow = Double.parseDouble(str2[i]);
			}
			if(codingPattern == CodingPattern.MSR){
				for(int i = 0; i < d-k+1; i++){
	//				System.out.println(str2[i]);
					V[al.get(i).getKey()].flow = Double.parseDouble(str2[i]);
				}
				for(int i = d-k+1; i < d; i++){
					//System.out.println(al.get(i).getKey());
					V[al.get(i).getKey()].flow = V[al.get(i-1).getKey()].flow;
				}
			}
		}catch(Exception e){
			e.printStackTrace();
		}
//		double bandwidthSum = 0;
//		for(int i = 0; i < d-k+1; i++)
//			bandwidthSum += al.get(i).getValue();
//		for(int i = 0; i < d-k+1; i++)
//			V[al.get(i).getKey()].flow[flowPattern] = alpha*(al.get(i).getValue()/bandwidthSum);
//		for(int i = d-k+1; i < d; i++)
//			V[al.get(i).getKey()].flow[flowPattern] = V[al.get(i-1).getKey()].flow[flowPattern];
	}
	
	private void clearLink(){
		for(int i = 0; i < nodeNum; i++){
			V[i].flow = 0;
		}
	}
	
	private double CompRouter4(int u, int v){
//		Node[] V2 = V.clone();

		Node[] V2 = new Node[nodeNum];
		for(int i = 0; i < nodeNum; i++)
			V2[i] = (Node)V[i].clone();
		
		double time = 0;
		V2[u].pi = v;
		V2[u].flow = beta;
		while(u != r){
			int x = V2[u].pi;
			if(V2[u].flow + V2[x].flow <= alpha)
				V2[x].flow = V2[u].flow + V2[x].flow;
			else
				V2[x].flow = alpha;
			if(time < V2[u].flow/w[u][x]){
				time = V2[u].flow/w[u][x];
//				if(time > minTime)
//					return minTime;
			}
			u = x;
		}
		//System.out.println("\n"+"V");
		//Print_node(V);
		//System.out.println("V2");
		//Print_node(V2);
		//System.out.println("time = " + time);
		return time;
	}
	
	private LinkedList<Integer> generalMaxSpanningTree_Prim(){
		LinkedList<Integer> A = new LinkedList<Integer>();
		double[] key = new double[nodeNum];
		for(int i = 0; i < nodeNum; i++){
			key[i] = 0;
			V[i].pi = -1;
		}
		key[r] = MAX_VALUE*10;
		LinkedList<Integer> Q = new LinkedList<Integer>();
		for(int i = 0; i < nodeNum; i++){
			Q.add(Integer.valueOf(i));
		}
		//System.out.println(Q);
		//Print_node();
		while(Q.size() > 0){
			//System.out.println(Q);
			int u = extractMax(key, Q);
			//System.out.println("u = "+ u);
			//Q.remove(Integer.valueOf(u));
			A.add(Integer.valueOf(u));
			for(int v = 0; v < nodeNum; v++){
				if(Q.contains(v) && w[v][u] > key[v]){
					V[v].pi = u;
					key[v] = w[v][u];
				}
			}
		}
		A.poll();
		return A;
	}
	
	private int extractMax(double[] key, LinkedList<Integer> Q){
		double maxKey = -1;
		int u = -1;
		Iterator<Integer> itr = Q.iterator();
		while(itr.hasNext()){
			int t = itr.next().intValue();
			if(key[t] > maxKey){
				maxKey = key[t];
				u = t;
			}
		}
		Q.remove(Integer.valueOf(u));
		return u;
	}
	
	private void JunAdjustRootLink(int x){
		int[][] linkGraph = getTreeLinkGraph_W();
		int rootLinkNum = 0;
		for(int i = 0; i < nodeNum; i++){
			if(linkGraph[r][i] != 0)
				rootLinkNum++;
		}
		if(rootLinkNum >= x)
			return;
		else{
			HashMap<Integer, Double> map = new HashMap<Integer, Double>();
			for(int i = 0; i < nodeNum; i++)
				if(linkGraph[r][i] == 0)
					map.put(Integer.valueOf(i), Double.valueOf(w[i][r]));
			ArrayList<Entry<Integer, Double> >  al = HashMapLargeToSmallSort(map);
			//System.out.println(al);
			for(int i = 0; i < x - rootLinkNum; i++){
				int j = al.get(i).getKey().intValue();
				V[j].pi = r;
			}
		}
	}
	
	public int[][] getTreeLinkGraph_W(){
    	int[][] linkGraph = new int[n][n];
    	int len = involvedNodeSet.length;
    	//System.out.println("involvedNodeSet:" + Arrays.toString(involvedNodeSet));
    	for(int i = 0; i < len; i++){
    		for(int j = 0; j < len; j++){
    			if(i == V[j].pi){
    				int i2 = involvedNodeSet[i];
    				int j2 = involvedNodeSet[j];
    				linkGraph[i2][j2] = 1;
    			}
    		}
    	}
    	return linkGraph;
    }
	
	private ArrayList<Entry<Integer, Double> > HashMapLargeToSmallSort(HashMap<Integer, Double> map){
		ArrayList<Entry<Integer, Double> > al = new ArrayList<Entry<Integer, Double> > (map.entrySet());
		Collections.sort(al, new Comparator<Map.Entry<Integer, Double> >(){
			public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2){
				if (o2.getValue().doubleValue() < o1.getValue().doubleValue())
					return -1;
				else if(o2.getValue().doubleValue() == o1.getValue().doubleValue())
					return 0;
				else
					return 1;
			}
		});
		return al;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Tree T = new Tree(5,4,3,100, CodingPattern.MBR, 0);
		T.testIndirectedLinkBandwidth3();
		
		Link starMinTimeLink = T.starMinTime();
		T.Print_bottleneck(starMinTimeLink);
		T.Print_Node();
		
		Link starMinTimeLinkChgFlowLink = T.starMinTimeChgFlow();
		T.Print_bottleneck(starMinTimeLinkChgFlowLink);
		T.Print_Node();
		
		Link singleMinTimeTreeHeuristicIncludeNodeOneByOneLink = T.singleMinTimeTreeHeuristicIncludeNodeOneByOne();
		T.Print_bottleneck(singleMinTimeTreeHeuristicIncludeNodeOneByOneLink);
		T.Print_Node();
		
		Link singleMinTimeTreeHeuristicIncludeNodeOneByOneChgFlowLink = T.singleMinTimeTreeHeuristicIncludeNodeOneByOneChgFlow();
		T.Print_bottleneck(singleMinTimeTreeHeuristicIncludeNodeOneByOneChgFlowLink);
		T.Print_Node();
		
		Link JunSingleMaxSpanningTreeLink =  T.JunSingleMaxSpanningTree();
		T.Print_bottleneck(JunSingleMaxSpanningTreeLink);
		T.Print_Node();
	}
}
