package edu.uta.futureye.tutorial;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import edu.uta.futureye.algebra.SchurComplementStokesSolver;
import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.application.DataReader;
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
import edu.uta.futureye.lib.element.FiniteElementType;
import edu.uta.futureye.lib.weakform.WeakFormNavierStokes;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjIndex;
import edu.uta.futureye.util.container.ObjList;


/**
 * Problem: Navier-Stokes, flow around a cylinder
 * 
 * Ref.
 * 1. Alexander N. BROOKS Streamline Upwind/Petrow-Galerkin Formulations for Convection Dominated
 *    Flows With Particular Emphasis on the Incompressible Navier-Stokes Equations
 *    
 * 2. M. Schafer and S. Turek Benchmark Computations of Laminar Flow Arond a Cylinder
 * @author liuyueming
 *
 */
public class NavierStokesCylinder {
	protected String file = null;
	protected String outputFolder = "tutorial/NavierStokesCylinder";
	protected Mesh mesh = null;
	protected Mesh meshOld = null;
	
	//Navier-Stokes Weak Form (For Picard Iteration)
	protected WeakFormNavierStokes weakForm = new WeakFormNavierStokes();
	//Assembler
	protected AssemblerVector assembler = null;
	//Dirichlet boundary condition
	protected VectorFunction diri = null;
	//Previous Velocity
	protected VectorFunction U = new SpaceVectorFunction(2);
	
	//delta t
	protected double dt = 0.02;
	//viscosity
	protected double mu = 0.001; 
	
	FiniteElementType fe = null;
	
	int maxNonlinearIter = 30;
	double nonlinearError = 1e-2;
	int maxTimeStep = 1000;
	
	
	public void init(int testCaseNo) {
		//Read a triangle mesh from an input file
		if(testCaseNo==1 || testCaseNo ==2) {
			//file = "benchmark_cylinder1"; //triangle mesh
			file = "benchmark_cylinder2"; //triangle mesh
		} else if(testCaseNo == 3) {
			file = "benchmark_cylinder3"; //rectangle mesh
			//file = "benchmark_cylinder4"; //rectangle mesh
			//file = "benchmark_cylinder5"; //rectangle mesh
		} else {
			throw new FutureyeException("testCaseNo should be 1,2,3");
		}
		MeshReader reader = new MeshReader(file+".grd");
		MeshReader reader2 = new MeshReader(file+".grd");
		mesh = reader.read2DMesh();
		meshOld = reader2.read2DMesh();
		mesh.nVertex = mesh.getNodeList().size();
		
		//Add nodes for quadratic element
		if(testCaseNo==1 || testCaseNo==2) {
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
		if(testCaseNo == 1)
			fe = new FEQuadraticV_LinearP();//Quadratic Velocity - Linear Pressure Element
		else if(testCaseNo == 2)
			fe = new FEQuadraticV_ConstantP();
		else if(testCaseNo == 3)
			fe = new FEBilinearV_ConstantP();
		fe.initDOFIndexGenerator(mesh);
		for(int i=1;i<=eList.size();i++) {
			fe.assignTo(eList.at(i));
//			eList.at(i).printDOFInfo();
		}

		//Boundary condition
		diri = new SpaceVectorFunction(3);
		diri.set(1, new AbstractFunction("x","y") {
					@Override
					public double value(Variable v) {
						double x = v.get("x");
						double y = v.get("y");
						
						double H  = 0.41;
						double Um = 1.5; //0.3
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
			
		weakForm.setParam(FC.c(0.001),U,FC.c(1.0/dt));
		
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

//		Matrix C = stiff.getBlock(3, 3);
//		for(int i=1;i<=C.getRowDim();i++)
//			C.set(i, i, 0.000000001);
		
		assembler.imposeDirichletCondition(diri);
		System.out.println("Assemble done!");
		
		SchurComplementStokesSolver solver = 
			new SchurComplementStokesSolver(stiff,load);
		//solver.setCGInit(0.5);
		//solver.debug = true;
		return solver.solve();
		
	}
	
	public void run(int restartStep, int testCaseNo) {
		init(testCaseNo);
		
		if(restartStep>0) {
			Vector vecU = DataReader.readVector(String.format("./%s/%s_uv_final_t%02d.dat",
					outputFolder,file,restartStep),3);
			Vector vecV = DataReader.readVector(String.format("./%s/%s_uv_final_t%02d.dat",
					outputFolder,file,restartStep),4);
			U.set(1, new Vector2Function(vecU));
			U.set(2, new Vector2Function(vecV));
		} else {
			U.set(1, FC.c0);
			U.set(2, FC.c0);
		}
		BlockVector u = null;
		SpaceVectorFunction uk = new SpaceVectorFunction(2);
		for(int time=restartStep+1;time<this.maxTimeStep;time++) {
			System.out.println(">>>>>>>>>>>>>>>>>>>time="+time);
			uk.set(1, U.get(1));
			uk.set(2, U.get(2));
			for(int iter=0;iter<this.maxNonlinearIter;iter++) {
				
				u = nonlinearIter(time, iter, uk);
				
				//Compute norm of delta_u (not delta_v)
				int dim = u.getBlock(1).getDim();
				SparseVector delta_u = new SparseVector(dim);
				for(int i=1;i<=dim;i++)
					delta_u.set(i, 
							u.getBlock(1).get(i)-
							U.get(1).value(new Variable().setIndex(i)));
				
				U.set(1, new Vector2Function(u.getBlock(1)));
				U.set(2, new Vector2Function(u.getBlock(2)));
	
				System.out.println("Iter="+iter+" Error Norm2 (||u1_k+1 - u1_k||) = "+delta_u.norm2());
				
				if(delta_u.norm2() < this.nonlinearError) {
					Tools.plotVector(mesh, outputFolder, String.format("%s_uv_final_t%03d.dat",file,time), 
							u.getBlock(1), u.getBlock(2));
					Tools.plotVector(meshOld, outputFolder, String.format("%s_p_final_t%03d.dat",file,time), 
							Tools.valueOnElement2Node(mesh,u.getBlock(3)));
					break;
				} else {
					Tools.plotVector(mesh, outputFolder, String.format("%s_uv%02d_%02d.dat",file,time,iter), 
							u.getBlock(1), u.getBlock(2));
					Tools.plotVector(meshOld, outputFolder, String.format("%s_p%02d_%02d.dat",file,time,iter), 
							Tools.valueOnElement2Node(mesh,u.getBlock(3)));
				}
			}
		}
	}
	
	public void computeKeyValues(int step,Vector vecU, Vector vecV, Vector vecP) {
		final double[] cylinderCenter = {0.2,0.2};
		double cylinderRadius = 0.05;
		//Node frontPoint  = new Node(0,0.15,0.2);
		//Node endPoint    = new Node(0,0.25,0.2);
		//Node frontPoint  = new Node(0, 0.145557922, 0.2000718387);  //GN547( 0.145557922 0.2000718387 )
		//Node endPoint    = new Node(0, 0.2544421267, 0.1999237695);  //GN377( 0.2544421267 0.1999237695 )
		Node frontPoint  = new Node(0, 0.1395858728, 0.2002194199);  //GN546( 0.1395858728 0.2002194199
		Node endPoint    = new Node(0, 0.2604143184, 0.1998722827);  //GN378( 0.2604143184 0.1998722827
		final SpaceVector v1 = new SpaceVector(2);
		final SpaceVector v2 = new SpaceVector(2);
		v1.set(1, frontPoint.coord(1)-cylinderCenter[0]);
		v1.set(2, frontPoint.coord(2)-cylinderCenter[1]);
		
		NodeList nodes = mesh.getNodeList();
		int N = nodes.size();
		NodeList nodesOnCylinder = new NodeList();
		for(int i=1;i<=N;i++) {
			Node node = nodes.at(i);
			double x = node.coord(1);
			double y = node.coord(2);
			v2.set(1, x-cylinderCenter[0]);
			v2.set(2, y-cylinderCenter[1]);
			//if(Math.abs(v2.norm2()-cylinderRadius) < Constant.meshEps) {
			//if(0.058-v2.norm2()>0 && v2.norm2() - 0.053>0) {
			if(0.067-v2.norm2()>0 && v2.norm2() - 0.058>0) {
				nodesOnCylinder.add(node);
				if(node.coordEquals(frontPoint)) {
					frontPoint.globalIndex = node.globalIndex;
				}
				if(node.coordEquals(endPoint)) {
					endPoint.globalIndex = node.globalIndex;
				}
			}
		}
		//排序结点从frontPoint开始，顺时针排序
		List<Node> l = nodesOnCylinder.toList();
		Collections.sort(l, new Comparator<Node>() {
			@Override
			public int compare(Node o1, Node o2) {
				v2.set(1, o1.coord(1)-cylinderCenter[0]);
				v2.set(2, o1.coord(2)-cylinderCenter[1]);
				double dot1 = v1.dot(v2);
				
				v2.set(1, o2.coord(1)-cylinderCenter[0]);
				v2.set(2, o2.coord(2)-cylinderCenter[1]);
				double dot2 = v1.dot(v2);
				
				if(o1.coord(2) >= cylinderCenter[1] && o2.coord(2) >= cylinderCenter[1]) {
					return (dot1-dot2)>0?-1:1;
				} else if(o1.coord(2) >= cylinderCenter[1] && o2.coord(2) < cylinderCenter[1]) {
					return -1;
				} else if(o1.coord(2) < cylinderCenter[1] && o2.coord(2) >= cylinderCenter[1]) {
					return 1;
				} else {
					return (dot1-dot2)>0?1:-1;
				}
			}
		});
//		for(int i=1;i<=nodesOnCylinder.size();i++) {
//			System.out.println(nodesOnCylinder.at(i));
//		}
		
		double H = 0.41; //the channel height
		double U = 1.5; //initial velocity at (0,H/2)
		double Ubar = 2.0*U/3.0; //the mean velocity
		double rho = 1.0; //fluid density
		double mu = 0.001; //kinematic viscosity
		double D = 0.1; //cylinder diameter
		double Re = Ubar*D/mu;
		
		double last_vtx=0.0, last_vty=0.0;
		Node last_node = null;
		double FD = 0.0;
		double FL = 0.0;
		for(int i=1;i<=nodesOnCylinder.size()+1;i++) {
			int index = i%nodesOnCylinder.size();
			if(index==0) index=nodesOnCylinder.size();
			Node node = nodesOnCylinder.at(index);
			v1.set(1, node.coord(1)-cylinderCenter[0]);
			v1.set(2, node.coord(2)-cylinderCenter[1]);
			
			//outer normal vector (nx ny)'
			double nx = v1.get(1)/v1.norm2();
			double ny = v1.get(2)/v1.norm2();
			//tangent vector (tx ty)'
			double tx = ny;
			double ty = -nx;
			//tangential velocity
			double u = vecU.get(node.globalIndex);
			double v = vecV.get(node.globalIndex);
			double vtLen = u*tx+v*ty;
			double vtx = vtLen*tx;
			double vty = vtLen*ty;
			
			if(i>=2) {
				double dSx = node.coord(1)-last_node.coord(1);
				double dSy = node.coord(2)-last_node.coord(2);
				double dS = Math.sqrt(dSx*dSx+dSy*dSy);
//				double dvtx = (vtx-last_vtx)/dSx;
//				double dvty = (vty-last_vty)/dSy;
				double dn=v1.norm2() - cylinderRadius;//
				double dvtx = vtx/dn;
				double dvty = vty/dn;		
				double P = vecP.get(node.globalIndex);
//				double vt_n = dvtx*nx+dvty*ny;
				double vt_n = Math.sqrt(dvtx*dvtx+dvty*dvty);
				FD += (rho*mu*vt_n*ny-P*nx)*dS;
				FL += -(rho*mu*vt_n*nx+P*ny)*dS;
				//System.out.println(FD+" "+FL);
			}
			
			last_vtx = vtx;
			last_vty = vty;
			last_node = node;
		}
		double cD = 2*FD/(rho*Ubar*Ubar*D);
		double cL = 2*FL/(rho*Ubar*Ubar*D);
		double f = 1.0/(9.0*0.05);//1/(step*dt)通过观察cD或cL来得出涡街频率
		double St = D*f/Ubar;
		double DeltaP = vecP.get(frontPoint.globalIndex)-vecP.get(endPoint.globalIndex);
		
		//System.out.println("index cD cL St DeltaP");
		System.out.println(String.format("%d\t%f\t%f\t%f\t%f", step,cD,cL,St,DeltaP));
		
	}
	
	public void computeKeyValues(int step) {
		Vector vecU = DataReader.readVector(String.format("./%s/%s_uv_final_t%02d.dat",
				outputFolder,file,step),3);
		Vector vecV = DataReader.readVector(String.format("./%s/%s_uv_final_t%02d.dat",
				outputFolder,file,step),4);
		Vector vecP = DataReader.readVector(String.format("./%s/%s_p_final_t%02d.dat",
				outputFolder,file,step));
		computeKeyValues(step,vecU,vecV,vecP);
	}
	
	public static void main(String[] args) {
		NavierStokesCylinder NS = new NavierStokesCylinder();
		NS.dt = 0.02;
		int timeStep = 0;
		int testCaseNo = 3;
		if(args.length == 1)
			NS.dt = Double.parseDouble(args[0]);
		if(args.length == 2)
			timeStep = Integer.parseInt(args[1]);
		if(args.length == 3)
			testCaseNo = Integer.parseInt(args[2]);

		System.out.println("dt="+NS.dt);
		System.out.println("timeStep="+timeStep);
		System.out.println("testCaseNo="+timeStep);
		
		System.out.println("mu="+NS.mu);
		System.out.println("maxTimeStep="+NS.maxTimeStep);
		System.out.println("maxNonlinearIter="+NS.maxNonlinearIter);
		System.out.println("nonlinearError="+NS.nonlinearError);
		
		NS.run(timeStep,testCaseNo);
		
//		//computeKeyValues
//		NS.init(3);
//		System.out.println("index cD cL St DeltaP");
//		for(int i=1;i<=620;i++)
//			NS.computeKeyValues(i);
	}
}