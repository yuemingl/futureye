package edu.uta.futureye.application;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MeshReader;

public class KlibanovNew {
	String outputFolder = "KlibanovNew";
	Mesh omega = null;
	Mesh omega0 = null;
	
	/**
	 * Load mesh
	 */
	public void init() {
		String gridFileBig = "KlibanovNewOmega0.grd";
		String gridFileSmall = "KlibanovNewOmega.grd";
		MeshReader readerForward = new MeshReader(gridFileBig);
		omega0 = readerForward.read2DMesh();
		MeshReader readerGCM = new MeshReader(gridFileSmall);
		omega = readerGCM.read2DMesh();
		
		//Use element library to assign degree of freedom (DOF) to element
		Tools.assignLinearShapFunction(omega0);
		Tools.assignLinearShapFunction(omega);
		omega0.computeNodeBelongsToElements();
		omega0.computeNeighborNodes();
		omega.computeNodeBelongsToElements();
		omega.computeNeighborNodes();

	}
	
	/**
	 * Laplace transform
	 * 
	 * U(x,z,s) = ReU + i*ImU = int_{-\inf}^{\inf}{u(x,z,x0)*e^{isx0}}dx0
	 * 
	 * where i=sqrt(-1)
	 * 
	 * @param u
	 * @param ReU
	 * @param ImU
	 */
	public static void transformLaplace(
			Vector u, double s,  //input
			Vector ReU, Vector ImU //output
			) {
		//write your code here
		
	}
	
	public Vector solveForwardProblem() {
		//Parameters of forward problem
		double k = 3.0;
		Function ax = ModelParam.getMu_a(2.0, 2.6, 0.5, //center: (x,y;r)
									6.0, //maxMu_a: 6.0
									1); //type: one inclusion
		
		ModelDOT model = new ModelDOT();
		model.setMu_sp(1.0/3.0);
		model.setMu_a(ax.A(k)); //k+a(x) 
		model.setLightPosition(-0.5, 4.0);
		Tools.plotFunction(omega0, outputFolder, "mu_a_real.dat", model.getMu_a());
		
		Vector u = model.solveNeumann(omega0);
		Tools.plotVector(omega0, outputFolder, "u.dat", u);
		return u;
	}
	
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		KlibanovNew prm = new KlibanovNew();
		prm.init();
		Vector u0 = prm.solveForwardProblem(); //forward solution on mesh omega0
		Vector u = Tools.extractData(prm.omega0, prm.omega, u0); //forward solution on mesh omega
		
		
		Vector ReU = new SpaceVector(u.getDim());
		Vector ImU = new SpaceVector(u.getDim());
		
		double s = 1;
		transformLaplace(u,s,ReU,ImU);
		
		//update tail
		//write you code here to update tail
		
		//solve q
		//Please modify the file GCMModelNew.java to implement 
		//the new algorithm based on the current implementation of GCM
		GCMModelNew model = new GCMModelNew(prm.outputFolder);
		//model.solveGCM(mesh, N, s, phi, tailT)
	}

}
