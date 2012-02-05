package edu.uta.futureye.application;

import java.util.HashMap;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.FaceLocal;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.geometry.GeoEntity3D;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.assembler.AssemblerScalarFast;
import edu.uta.futureye.lib.element.FETrilinearHexahedron;
import edu.uta.futureye.lib.weakform.WeakFormLaplace;
import edu.uta.futureye.lib.weakform.WeakFormLaplace3D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjList;
import static edu.uta.futureye.function.operator.FMath.*;

public class HumanRealSimulate {
	String workFolder = "HumanReal";
	//light source position
	public Variable lightPosition = null; 
	Function delta = null;
	Mesh meshOmega = null;
	
	double mu_sp = 0.0;
	Function mu_a = null;
	
	public void setModelParam(double mu_sp, Function mu_a) {
		this.mu_sp = mu_sp;
		this.mu_a = mu_a;
	}
	public void setLightPosition(double x,double y,double z) {
		this.lightPosition = new Variable();
		this.lightPosition.set("x", x);
		this.lightPosition.set("y", y);
		this.lightPosition.set("z", z);
		delta = new FDelta(this.lightPosition,0.01,2e5);
	}
	
	/**
	 * Solver the following model problem:
	 * 
	 *   -\nabla{(1/k)*\nabla{u}} + a*u = f
	 * 
	 * where
	 *   f = delta(x-x0)*delta(y-y0)*delta(z-z0)
	 *   u = u(x,y,z)
	 *   k = 1/(3*mu_s')
	 *   a = a(x,y,z) = mu_a(x,y,z)
	 */
	public Vector solveForward(String meshName, String info) {
		System.out.println("Solve forward: mesh \""+meshName + "\" info:" + info);
		
		MeshReader reader = new MeshReader(workFolder+"/"+meshName+".grd");
		meshOmega = reader.read3DMesh(); //3D
		meshOmega.computeNodeBelongsToElements(); //worked in 3D
		meshOmega.computeGlobalEdge();
		
		HashMap<NodeType, Function> mapNTF = new HashMap<NodeType, Function>();
		mapNTF.put(NodeType.Neumann, null);
		meshOmega.markBorderNode(mapNTF);
		meshOmega.writeNodesInfo(workFolder+"/meshInfo.dat");
		
        ElementList eList = meshOmega.getElementList();
        FETrilinearHexahedron feTLH = new FETrilinearHexahedron();
        for(int i=1;i<=eList.size();i++)
        	feTLH.assignTo(eList.at(i));

//		Tools.plotFunction(meshOmega, workFolder, "mu_a.dat", this.mu_a);
//		Tools.plotFunction(meshOmega, workFolder, "delta.dat", this.delta);
		
//		int count=0,count2=0;
//		for(Element e : mesh.getElementList()) {
//			if(e.isBorderElement()) count++;
//			GeoEntity3D entity = e.getGeoEntity3D();
//				ObjList<FaceLocal> faces = entity.getFaces();
//				for(int i=1;i<=faces.size();i++) {
//					if(faces.at(i).isBorderFace())
//						count2++;
//				}
//			//检查边界单元的结点是否在同一个面上：有一个坐标值全相同（因为网格是正规的）
//			ElementList beList = e.getBorderElements();
//			for(Element be : beList) {
//				boolean find2 = false;
//				for(int j=1;j<=3;j++) {
//					boolean find = true;
//					for(int i=2;i<=be.nodes.size();i++) {
//						if(be.nodes.at(1).coord(j)!=be.nodes.at(i).coord(j)) {
//							find = false;
//							break;
//						}
//					}
//					if(find) { find2 = true; break; }
//				}
//				if(!find2) {
//					System.out.println("Error : "+be);
//				}
//			}
//		}
//		System.out.println(count);
//		System.out.println(count2);
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace3D weakForm = new WeakFormLaplace3D();
		
		weakForm.setF(delta);
		//weakForm.setF(C(1));
		
		weakForm.setParam(
				C(1./(3*this.mu_sp)), 
				this.mu_a, 
				null, 
				C(1./(3*this.mu_sp))//null
			);
		
		Assembler assembler = new AssemblerScalar(meshOmega, weakForm);
		System.out.println("Begin Assemble...");
		long begin = System.currentTimeMillis();
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		System.out.println("Assemble done!");
		long end = System.currentTimeMillis();
		System.out.println("Time used:"+(end-begin));
		System.out.println("imposeDirichletCondition()...");
		begin = System.currentTimeMillis();
		assembler.imposeDirichletCondition(FC.C0);
		end = System.currentTimeMillis();
		System.out.println("imposeDirichletCondition() time used:"+(end-begin));
		
		SolverJBLAS solver = new SolverJBLAS();
		begin = System.currentTimeMillis();
		Vector u = solver.solveDGESV(stiff, load);
		end = System.currentTimeMillis();
		System.out.println("SolverJBLAS time used:"+(end-begin));
	    
	    //Tools.plotVector(meshOmega, workFolder, meshName+"_u.dat", u);
		return u;
	}
	
	/**
	 * Measure surface is at y=0 of Omega
	 * Measure surface(xi,yi) <==> Omega(yi,0,xi)
	 * 
	 * @param measureSurface
	 * @param u
	 * @return
	 */
	public Vector extractMeasureSurfaceValues(Mesh measureSurface, Vector u) {
		NodeList surfNodes = measureSurface.getNodeList();
		Vector rlt = new SpaceVector(surfNodes.size());
		double []coords = new double[3];
		for(int i=1;i<=surfNodes.size();i++) {
			Node sn = surfNodes.at(i);
			coords[0] = sn.coord(2);
			coords[1] = 0.0;
			coords[2] = sn.coord(1);
			Node on = this.meshOmega.containNode(new Node().set(0,coords));
			if(on != null) 
				rlt.set(i, u.get(on.globalIndex));
		}
		return rlt;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		HumanRealSimulate sim = new HumanRealSimulate();
		sim.setLightPosition(1.5,-0.5,4);
		//sim.setLightPosition(0.5,-0.5,15);
		//sim.setLightPosition(0.5,-0.5,17.5);
		
		double mu_sp = 10;
		Function fmu_a = new AbstractFunction("x","y","z"){
			double x0 = 5,y0=-1.5, z0=4;//z0=13
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				double z = v.get("z");
				if(Math.sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))<1.0)
					return 0.5;
				else
					return 0.1;
			}
		};
		Function fmu_a0 = C(0.1);
		
		String meshName = "mesh3DLeft";
		//String meshName = "mesh3DRight";
		String meshMeasureSurfaceName = "meshMeasureSurfaceLeft";
		//String meshMeasureSurfaceName = "meshMeasureSurfaceRight";
		MeshReader reader = new MeshReader(sim.workFolder+"/"+meshMeasureSurfaceName+".grd");
		Mesh meshMeasureSurface = reader.read2DMesh(); //2D
		
		sim.setModelParam(mu_sp, fmu_a0);
		Vector u0 = sim.solveForward(meshName,"background");
		Tools.plotVector(sim.meshOmega, sim.workFolder, meshName+"_u0.dat", u0);
		Vector meaSurf0 = sim.extractMeasureSurfaceValues(meshMeasureSurface, u0);
		Tools.plotVector(meshMeasureSurface, sim.workFolder, meshName+"_u0_MeaSurf.dat", meaSurf0);

		sim.setModelParam(mu_sp, fmu_a);
		Vector u = sim.solveForward(meshName,"inclusion");
		Tools.plotVector(sim.meshOmega, sim.workFolder, meshName+"_u.dat", u);
		Vector meaSurf = sim.extractMeasureSurfaceValues(meshMeasureSurface, u);
		Tools.plotVector(meshMeasureSurface, sim.workFolder, meshName+"_u_MeaSurf.dat", meaSurf);
		
		double []setup_x = {1,3,5,7,  2,  4,  6,1,3,5,7};
		double []setup_y = {3,3,3,3,4.5,4.5,4.5,6,6,6,6};
		//meaSurf - meaSurf0
		Vector2Function fmeaSurf0 = new Vector2Function(meaSurf0, meshMeasureSurface, "x","y");
		Vector2Function fmeaSurf = new Vector2Function(meaSurf, meshMeasureSurface, "x","y");
		Vector diffIntensity = new SpaceVector(setup_x.length);
		Vector diffOD = new SpaceVector(setup_x.length);
		for(int i=0;i<setup_x.length;i++) {
			Variable v = new Variable("x",setup_x[i]).set("y",setup_y[i]);
			diffIntensity.set(i+1, fmeaSurf.value(v)-fmeaSurf0.value(v));
			diffOD.set(i+1, -Math.log(fmeaSurf0.value(v)/fmeaSurf.value(v)));
		}
		
		MeshWriter.writeTechplotLine(sim.workFolder + "/detector.txt", 
				new SpaceVector(setup_x), 
				new SpaceVector(setup_y),
				diffIntensity,
				diffOD
				);
		
	}

}
