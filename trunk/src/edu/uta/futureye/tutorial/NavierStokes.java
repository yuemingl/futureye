package edu.uta.futureye.tutorial;

import java.util.HashMap;

import edu.uta.futureye.algebra.SchurComplementStokesSolver;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.Vertex;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.lib.assembler.AssemblerVector;
import edu.uta.futureye.lib.element.FEBilinearV_ConstantP;
import edu.uta.futureye.lib.element.FEQuadraticV_ConstantP;
import edu.uta.futureye.lib.element.FEQuadraticV_LinearP;
import edu.uta.futureye.lib.weakform.WeakFormNavierStokes;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;


/**
 * Problem: Navier-Stokes
 * 
 * @author liuyueming
 *
 */
public class NavierStokes {
	//protected String file = "benchmark_cylinder1"; //triangle mesh
	//protected String file = "benchmark_cylinder2"; //triangle mesh
	protected String file = "benchmark_cylinder3"; //rectangle mesh
	protected static String outputFolder = "tutorial/NavierStokes";
	protected Mesh mesh = null;
	protected Mesh meshOld = null;
	
	//Quadratic Velocity - Linear Pressure Element
	int feFlag = 3;
	//protected FEQuadraticV_LinearP fe = new FEQuadraticV_LinearP(); //feFlag=1
	//protected FEQuadraticV_ConstantP fe = new FEQuadraticV_ConstantP(); //feFlag=2
	protected FEBilinearV_ConstantP fe = new FEBilinearV_ConstantP(); //feFla=3
	
	//Navier-Stokes Weak Form (For Picard Iteration)
	protected WeakFormNavierStokes weakForm = new WeakFormNavierStokes();
	//Assembler
	protected AssemblerVector assembler = null;
	//Dirichlet boundary condition
	protected VectorFunction diri = null;
	//Previous Velocity
	protected VectorFunction U = new SpaceVectorFunction(2);
	
	//delta t
	protected double dt = 0.05;
	
	public void init() {
		//Read a triangle mesh from an input file
		MeshReader reader = new MeshReader(file+".grd");
		MeshReader reader2 = new MeshReader(file+".grd");
		mesh = reader.read2DMesh();
		meshOld = reader2.read2DMesh();
		mesh.nVertex = mesh.getNodeList().size();
		
		//Add nodes for quadratic element
		if(feFlag==1 && feFlag==2) {
			for(int i=1;i<=mesh.getElementList().size();i++) {
				Element e = mesh.getElementList().at(i);
				e.adjustVerticeToCounterClockwise();
				ObjList<EdgeLocal> edges = e.edges();
				int nNode = e.nodes.size();
				for(int j=1;j<=edges.size();j++) {
					EdgeLocal edge = edges.at(j);
					Vertex l = edge.beginVertex();
					Vertex r = edge.endVertex();
					double cx = (l.coord(1)+r.coord(1))/2.0;
					double cy = (l.coord(2)+r.coord(2))/2.0;
					Node node = new Node(mesh.getNodeList().size()+1, cx,cy);
					Node findNode = mesh.containNode(node);
					if(findNode == null) {
						edge.addEdgeNode(new NodeLocal(++nNode,node));
						mesh.addNode(node);
					} else {
						edge.addEdgeNode(new NodeLocal(++nNode,findNode));
					}
				}
				e.applyChange();
			}
		}
		
		//Geometry relationship
		mesh.computeNodeBelongsToElements();
		
		ElementList eList = mesh.getElementList();
		NodeList nodes = mesh.getNodeList();
//		for(int i=1;i<=eList.size();i++) {
//			System.out.println(i+"  " + eList.at(i));
//		}
		
		//Mark border type
		//指定u,v边界
		HashMap<NodeType, Function> mapNTF_uv = new HashMap<NodeType, Function>();
		mapNTF_uv.put(NodeType.Dirichlet, new AbstractFunction("x","y") {
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				//upside and down side
				if(y < Constant.meshEps || Math.abs(y-0.41) < Constant.meshEps)
					return 1;
				//left side
				else if(x < Constant.meshEps)
					return 1;
				
				//cylinder
				if(Math.sqrt((x-0.2)*(x-0.2)+(y-0.2)*(y-0.2)) <= 0.05+Constant.meshEps)
					return 1;
				
				return 0;
			}
		});
		//u,v其他边界
		mapNTF_uv.put(NodeType.Neumann, null);
		
		//指定p边界
		HashMap<NodeType, Function> mapNTF_p = new HashMap<NodeType, Function>();
		mapNTF_p.put(NodeType.Dirichlet, new AbstractFunction("x","y") {
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				//right side p=0
				if(Math.abs(x-2.2) < Constant.meshEps)
					return 1;
				else
					return 0;
			}
		});
		//p其他边界
		mapNTF_p.put(NodeType.Neumann, null);
		
		mesh.markBorderNode(new ObjIndex(1,2),mapNTF_uv);
		mesh.markBorderNode(3,mapNTF_p);

		//Use element library to assign degree of freedom (DOF) to element
		for(int i=1;i<=eList.size();i++) {
			System.out.println(i+"  " + eList.at(i));
		}

		fe.initDOFIndexCounter(nodes.size());
		for(int i=1;i<=eList.size();i++) {
			fe.assignTo(eList.at(i));
			eList.at(i).printDOFInfo();
		}

		//Boundary condition
		diri = new SpaceVectorFunction(3);
		diri.set(1, new AbstractFunction("x","y") {
					@Override
					public double value(Variable v) {
						double x = v.get("x");
						double y = v.get("y");
						
						double H  = 0.41;
						double Um = 0.3;
						//left boundary
						if(x < Constant.meshEps)
							return 4*Um*y*(H-y)/(H*H);
						else
							return 0.0;
					}
				});
		diri.set(2, FC.c0);
		diri.set(3, FC.c0);
	}
	
	public BlockVector nonlinearIter(int time, int nIter, SpaceVectorFunction uk) {
		//Right hand side(RHS): f = (0,0)'
		if(time==0)
			weakForm.setF(new SpaceVectorFunction(FC.c0,FC.c0));
		else
			weakForm.setF(new SpaceVectorFunction(uk.get(1).D(dt),uk.get(2).D(dt)));
			
		weakForm.setParam(FC.c(0.005),U,FC.c(1.0/dt));
		
		//Robin:  k*u_n + d*u - p\mathbf{n} = 0
		VectorFunction d = new SpaceVectorFunction(2);
		d.set(1, FC.c0);
		d.set(2, FC.c0);
		weakForm.setRobin(d);
		
		assembler = new AssemblerVector(mesh, weakForm,fe);
		System.out.println("Begin Assemble...");
		assembler.assemble();
		BlockMatrix stiff = (BlockMatrix)assembler.getStiffnessMatrix();
		BlockVector load = (BlockVector)assembler.getLoadVector();
		//stiff.print();
		//load.print();
		//System.out.println(load.norm2());

		Matrix C = stiff.getBlock(3, 3);
		for(int i=1;i<=C.getRowDim();i++)
			C.set(i, i, 0.00001);
		
		assembler.imposeDirichletCondition(diri);
		//load.getBlock(1).print();
		//load.getBlock(2).print();
		//load.getBlock(3).print();
		System.out.println("Assemble done!");
		
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		
		return solver.solve();
		
	}
	
	public void run() {
		init();
		
		U.set(1, FC.c(0.1));
		U.set(2, FC.c0);
		BlockVector u = null;
		SpaceVectorFunction uk = new SpaceVectorFunction(2);
		for(int time=0;time<200;time++) {
			System.out.println(">>>>>>>>>>>>>>>>>>>time="+time);
			uk.set(1, U.get(1));
			uk.set(2, U.get(2));
			for(int iter=0;iter<30;iter++) {
				
				u = nonlinearIter(time, iter, uk);
				
				int dim = u.getBlock(1).getDim();
				SparseVector tmp = new SparseVector(dim);
				for(int i=1;i<=dim;i++)
					tmp.set(i, 
							u.getBlock(1).get(i)-
							U.get(1).value(new Variable().setIndex(i)));
				
				U.set(1, new Vector2Function(u.getBlock(1)));
				U.set(2, new Vector2Function(u.getBlock(2)));
	
//				System.out.println("u=");
//				for(int i=1;i<=u.getDim();i++)
//					System.out.println(String.format("%.3f", u.get(i)));
				System.out.println("Iter="+iter+" Error Norm2 (||u1_k+1 - u1_k||) = "+tmp.norm2());
				if(tmp.norm2() < 1e-5) {
					Tools.plotVector(mesh, outputFolder, String.format("%s_uv_final_t%02d.dat",file,time), 
							u.getBlock(1), u.getBlock(2));
					Tools.plotVector(meshOld, outputFolder, String.format("%s_p_final_t%02d.dat",file,time), 
							element2node(u.getBlock(3)));
					break;
				} else {
					Tools.plotVector(mesh, outputFolder, String.format("%s_uv%02d_%02d.dat",file,time,iter), 
							u.getBlock(1), u.getBlock(2));
					Tools.plotVector(meshOld, outputFolder, String.format("%s_p%02d_%02d.dat",file,time,iter), 
							element2node(u.getBlock(3)));
				}
			}
		}
	}
	
	//Convert values of pressure p from piecewise constant on element to node if necessary 
	public Vector element2node(Vector p) {
		if(p.getDim()==mesh.getElementList().size() && p.getDim() != mesh.getNodeList().size()) {
		    Vector pOnElement = p;
		    ElementList eList = mesh.getElementList();
		    Vector pOnNode = new SparseVector(mesh.getNodeList().size());
		    for(int i=1;i<=eList.size();i++) {
		    	NodeList nList = eList.at(i).nodes;
		    	for(int j=1;j<=nList.size();j++) {
		    		if(pOnElement.get(eList.at(i).globalIndex) > pOnNode.get(nList.at(j).globalIndex))
		    		pOnNode.set(
		    				nList.at(j).globalIndex,
		    				pOnElement.get(eList.at(i).globalIndex)
		    				);
		    	}
		    }
		    return pOnNode;
		} else {
			return p;
		}
	}
	
	public static void main(String[] args) {
		NavierStokes NS = new NavierStokes();
		NS.run();
	}
}
