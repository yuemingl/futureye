package edu.uta.futureye.core;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.util.EdgeList;
import edu.uta.futureye.util.ElementList;
import edu.uta.futureye.util.FutureEyeException;
import edu.uta.futureye.util.NodeList;

public class Refiner {

	/**
	 * ����Ƿ���Ҫ��ʽ������Χ���񣬵����еڶ��μ��ܵ�Ԫʱ��Ϊ�˱���һ���߽���ֻ��һ��hanging node��
	 * ������Ҫ���ܵ�ǰ��Ԫ��Χ�ĵ�Ԫ
	 * @param mesh
	 * @param eToRefine
	 */
	public static ElementList checkNeighborRefinement(ElementList eToRefine) {
		Map<Integer, Element> map = new HashMap<Integer,Element>();
		ElementList findNeighbors = new ElementList();
		
		for(int i=1;i<=eToRefine.size();i++) {
			Element e = eToRefine.at(i);
			ElementList eNeighbors = e.neighbors;
			for(int j=1;j<=eNeighbors.size();j++) {
				Element findNeighbor = eNeighbors.at(j);
				if(e.getLevel() - findNeighbor.getLevel()>=1) {
					map.put(findNeighbor.globalIndex, findNeighbor);
				}
			}
		}
		for(Element e : map.values())
			findNeighbors.add(e);
		return findNeighbors;
	}
	
	public static void directRefine(Mesh mesh, ElementList eToRefine) {
		ElementList meList = mesh.getElementList();
		NodeList mnList = mesh.getNodeList();
		
		for(int i=1;i<=eToRefine.size();i++) {
			Element e = eToRefine.at(i);
			//��ֹ�ظ�refine
			if(e.isRefined())
				continue;
			
			ElementList eList = new ElementList();
			List<Vertex> vList = e.getVertexList();
			
			if(vList.size() == 3) {
				/*
				 * 3
				 * | \
				 * |  \
				 * 6---5
				 * |\  |\
				 * | \ | \
				 * 1--4---2
				 */
				Vertex v1 = vList.get(0);
				Vertex v2 = vList.get(1);
				Vertex v3 = vList.get(2);
				
				Node node4 = new NodeRefined(2);
				node4.setCoord(1, (v1.coord(1)+v2.coord(1))/2.0);
				node4.setCoord(2, (v1.coord(2)+v2.coord(2))/2.0);
				node4.setLevel(e.getLevel()+1);
				Node tmpNode = mesh.containNode(node4);
				if(tmpNode != null) node4 = tmpNode;
				
				Node node5 = new NodeRefined(2);
				node5.setCoord(1, (v2.coord(1)+v3.coord(1))/2.0);
				node5.setCoord(2, (v2.coord(2)+v3.coord(2))/2.0);
				node5.setLevel(e.getLevel()+1);
				tmpNode = mesh.containNode(node5);
				if(tmpNode != null) node5 = tmpNode;
				
				Node node6 = new NodeRefined(2);
				node6.setCoord(1, (v3.coord(1)+v1.coord(1))/2.0);
				node6.setCoord(2, (v3.coord(2)+v1.coord(2))/2.0);
				node6.setLevel(e.getLevel()+1);
				tmpNode = mesh.containNode(node6);
				if(tmpNode != null) node6 = tmpNode;
				
				Node node1 = e.nodes.at(v1.localIndex);
				Node node2 = e.nodes.at(v2.localIndex);
				Node node3 = e.nodes.at(v3.localIndex);
				
				Element e1 = new Element();
				Element e2 = new Element();
				Element e3 = new Element();
				Element e4 = new Element();
				
				e1.addNode(node1, true);
				e1.addNode(node4, true);
				e1.addNode(node6, true);
				e1.parent = e;
				e1.setLevel(e.getLevel()+1);
				
				e2.addNode(node4, true);
				e2.addNode(node2, true);
				e2.addNode(node5, true);
				e2.parent = e;
				e2.setLevel(e.getLevel()+1);

				e3.addNode(node5, true);
				e3.addNode(node3, true);
				e3.addNode(node6, true);
				e3.parent = e;
				e3.setLevel(e.getLevel()+1);

				e4.addNode(node4, true);
				e4.addNode(node5, true);
				e4.addNode(node6, true);
				e4.parent = e;
				e4.setLevel(e.getLevel()+1);

				eList.add(e1);
				eList.add(e2);
				eList.add(e3);
				eList.add(e4);
			} else if(vList.size() == 4) {
				/*
				 * 	4--7--3
				 *  |  |  |
				 *  8--9--6
				 *  |  |  |
				 *  1--5--2
				 */
				Vertex v1 = vList.get(0);
				Vertex v2 = vList.get(1);
				Vertex v3 = vList.get(2);
				Vertex v4 = vList.get(3);
				
				Node node5 = new NodeRefined(2);
				node5.setCoord(1, (v1.coord(1)+v2.coord(1))/2.0);
				node5.setCoord(2, (v1.coord(2)+v2.coord(2))/2.0);
				node5.setLevel(e.getLevel()+1);
				Node tmpNode = mesh.containNode(node5);
				if(tmpNode != null) node5 = tmpNode;
				
				Node node6 = new NodeRefined(2);
				node6.setCoord(1, (v2.coord(1)+v3.coord(1))/2.0);
				node6.setCoord(2, (v2.coord(2)+v3.coord(2))/2.0);
				node6.setLevel(e.getLevel()+1);
				tmpNode = mesh.containNode(node6);
				if(tmpNode != null) node6 = tmpNode;
				
				Node node7 = new NodeRefined(2);
				node7.setCoord(1, (v3.coord(1)+v4.coord(1))/2.0);
				node7.setCoord(2, (v3.coord(2)+v4.coord(2))/2.0);
				node7.setLevel(e.getLevel()+1);
				tmpNode = mesh.containNode(node7);
				if(tmpNode != null) node7 = tmpNode;
				
				Node node8 = new NodeRefined(2);
				node8.setCoord(1, (v4.coord(1)+v1.coord(1))/2.0);
				node8.setCoord(2, (v4.coord(2)+v1.coord(2))/2.0);
				node8.setLevel(e.getLevel()+1);
				tmpNode = mesh.containNode(node8);
				if(tmpNode != null) node8 = tmpNode;
				
				Node node9 = new NodeRefined(2);
				node9.setCoord(1, (node5.coord(1)+node7.coord(1))/2.0);
				node9.setCoord(2, (node5.coord(2)+node7.coord(2))/2.0);
				node9.setLevel(e.getLevel()+1);
				tmpNode = mesh.containNode(node9);
				if(tmpNode != null) node9 = tmpNode;
				
				Node node1 = e.nodes.at(v1.localIndex);
				Node node2 = e.nodes.at(v2.localIndex);
				Node node3 = e.nodes.at(v3.localIndex);
				Node node4 = e.nodes.at(v4.localIndex);
				
				Element e1 = new Element();
				Element e2 = new Element();
				Element e3 = new Element();
				Element e4 = new Element();
				
				e1.addNode(node1, true);
				e1.addNode(node5, true);
				e1.addNode(node9, true);
				e1.addNode(node8, true);
				e1.parent = e;
				e1.setLevel(e.getLevel()+1);
				
				e2.addNode(node5, true);
				e2.addNode(node2, true);
				e2.addNode(node6, true);
				e2.addNode(node9, true);
				e2.parent = e;
				e2.setLevel(e.getLevel()+1);

				e3.addNode(node9, true);
				e3.addNode(node6, true);
				e3.addNode(node3, true);
				e3.addNode(node7, true);
				e3.parent = e;
				e3.setLevel(e.getLevel()+1);

				e4.addNode(node8, true);
				e4.addNode(node9, true);
				e4.addNode(node7, true);
				e4.addNode(node4, true);
				e4.parent = e;
				e4.setLevel(e.getLevel()+1);

				eList.add(e1);
				eList.add(e2);
				eList.add(e3);
				eList.add(e4);
			} else {
				FutureEyeException ex = new FutureEyeException("Unsupported element type for refinement!");
				ex.printStackTrace();
			}
			
			for(int j=1;j<=eList.size();j++) {
				Element eNew = eList.at(j);
				if(eNew.globalIndex == 0) {
					eNew.globalIndex = meList.size()+1;
					meList.add(eNew);
					for(int k=1;k<=eNew.nodes.size();k++) {
						Node nNew = eNew.nodes.at(k);
						if(nNew.globalIndex == 0) {
							nNew.globalIndex = mnList.size()+1;
							mnList.add(nNew);
						}
					}
				}	
			}
			
			e.chlids = eList;
		}

	}

	
	public static void refineOnce(Mesh mesh, ElementList eToRefine) {
		ElementList eList = mesh.getElementList();
		while(true) {
			ElementList eNeighbors = checkNeighborRefinement(eToRefine);
			if(eNeighbors.size()>0) {
				directRefine(mesh,eNeighbors);
				computeHangingNode(eNeighbors);
				for(int iToRe=1;iToRe<=eNeighbors.size();iToRe++) {
					eList.remove(eNeighbors.at(iToRe));
				}
				mesh.computeNodesBelongToElement();
				mesh.computeNeiborNode();
				mesh.computeNeighborElement();
			} else {
				directRefine(mesh,eToRefine);
				computeHangingNode(eToRefine);
				for(int iToRe=1;iToRe<=eToRefine.size();iToRe++) {
					eList.remove(eToRefine.at(iToRe));
				}
				mesh.computeNodesBelongToElement();
				mesh.computeNeiborNode();
				mesh.computeNeighborElement();
				break;
			}
		}
	}
	
	//����hanging node
	public static void computeHangingNode(ElementList eToRefine) {
		for(int iToRe=1;iToRe<=eToRefine.size();iToRe++) {
			ElementList eChilds = eToRefine.at(iToRe).chlids;
			for(int i=1;i<=eChilds.size();i++) {
				Element eChild = eChilds.at(i);
				for(int j=1;j<=eChild.nodes.size();j++) {
					Node nNew = eChild.nodes.at(j);
					//!!! Node��level�뵥Ԫlevel��ͬʱ�����hanging node����ֵ������׼�����¼���
					if(nNew.getLevel() == eChild.getLevel()) {
						NodeRefined nRefined = (NodeRefined)nNew;
						nRefined.clearConstrainNode();
					}
				}
			}
		}
		
		for(int iToRe=1;iToRe<=eToRefine.size();iToRe++) {
			ElementList eChilds = eToRefine.at(iToRe).chlids;
			ElementList eNeighbors = eToRefine.at(iToRe).neighbors;
			//ѭ����Ԫ�е�ÿ���ӵ�Ԫ
			for(int i=1;i<=eChilds.size();i++) {
				Element eChild = eChilds.at(i);
				//ѭ��ÿ���ӵ�Ԫ�Ľ��
				for(int j=1;j<=eChild.nodes.size();j++) {
					Node nNew = eChild.nodes.at(j);
					if(nNew.getLevel() == eChild.getLevel()) {
						NodeRefined nRefined = (NodeRefined)nNew;
						//ѭ����Ԫ�����ڵ�Ԫ���жϸýڵ��Ƿ�Hanging node��
						//����Ǳ߽絥Ԫ�ı߽��ϼ��ܣ�nRefined��û��ConstrainNode��Ĭ��Ϊ��Hanging node
						for(int k=1;k<=eNeighbors.size();k++) {
							Element eNeighbor = eNeighbors.at(k);
							//!!! ���ڵ�Ԫû�б��Ϊ����  ���� ���ڵ�Ԫ�Ĳ��Ҫ���ڸõ�Ԫ
							if(!eNeighbor.isRefined() && eNeighbor.getLevel()<eChild.getLevel()) {
								EdgeList edges = eNeighbors.at(k).getEdgeList();
								for(int kk=1;kk<=edges.size();kk++) {
									//�жϽ���Ƿ��ڴ�Ԫ����
									if(edges.at(kk).isCoordOnEdge(nRefined.coords())) {
										NodeList endNodes = edges.at(kk).getEndNodes();
										//System.out.println("Hanging node:"+nRefined.globalIndex+
										//		" on "+endNodes.at(1).globalIndex+" "+endNodes.at(2).globalIndex);
										nRefined.addConstrainNode(endNodes.at(1));
										nRefined.addConstrainNode(endNodes.at(2));
									}
								}
							}
						}
					}
				}
			}
		}
	}

}