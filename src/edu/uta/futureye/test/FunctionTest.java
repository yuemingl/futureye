package edu.uta.futureye.test;

import java.util.ArrayList;
import java.util.List;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.application.Tools;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.AbstractVectorFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.basic.DuDx;
import edu.uta.futureye.function.basic.FAx;
import edu.uta.futureye.function.basic.FAxpb;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FLinear1D;
import edu.uta.futureye.function.basic.FPolynomial1D;
import edu.uta.futureye.function.basic.FX;
import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.VectorFunction;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import static edu.uta.futureye.function.operator.FMath.*;

public class FunctionTest {
	
	public static void constantTest() {
		Function c0 = FC.C0;
		Function c1 = FC.C1;
		Function fx = FX.fx;
		
		System.out.println(c0.A(c0));
		System.out.println(c1.A(c1));
		System.out.println(c1.A(c0));
		System.out.println(c0.A(c1));
		System.out.println();
		System.out.println(fx.A(c0));
		System.out.println(fx.A(c1));
		System.out.println(c0.A(fx));
		System.out.println(c1.A(fx));
		System.out.println(fx.A(fx));
		System.out.println();
		
		System.out.println(c0.S(c0));
		System.out.println(c1.S(c1));
		System.out.println(c1.S(c0));
		System.out.println(c0.S(c1));
		System.out.println();
		System.out.println(fx.S(c0));
		System.out.println(fx.S(c1));
		System.out.println(c0.S(fx));
		System.out.println(c1.S(fx));
		System.out.println(fx.S(fx));
		System.out.println();
		
		System.out.println(c0.M(c0));
		System.out.println(c1.M(c1));
		System.out.println(c1.M(c0));
		System.out.println(c0.M(c1));
		System.out.println();
		System.out.println(fx.M(c0));
		System.out.println(fx.M(c1));
		System.out.println(c0.M(fx));
		System.out.println(c1.M(fx));
		System.out.println(fx.M(fx));
		System.out.println();
		
		System.out.println(c0.D(c0));
		System.out.println(c1.D(c1));
		System.out.println(c1.D(c0));
		System.out.println(c0.D(c1));
		System.out.println();
		System.out.println(fx.D(c0));
		System.out.println(fx.D(c1));
		System.out.println(c0.D(fx));
		System.out.println(c1.D(fx));
		System.out.println(fx.D(fx));
		System.out.println();
		
	}
	
	public static void test() {
		
		FAx fax1 = new FAx(3.0);
		System.out.println("fax1 = "+fax1);
		FAx fax2 = new FAx(1.0);
		System.out.println("fax2 = "+fax2);
		FAx fax3 = new FAx(0.0);
		System.out.println("fax3 = "+fax3);
		
		FAxpb faxpb1 = new FAxpb(3.0,2.0);
		System.out.println("faxpb1 = "+faxpb1);
		
		FC c0 = new FC(0.0);
		FC c1 = new FC(1.0);
		FC c2 = new FC(2.0);
		FC c3 = new FC(3.0);
		
		FX fi = FX.fx;
		Function ftest1 = c3.M(fi).A(c2);
		System.out.println("ftest1 = "+ftest1);
		Function ftest2 = c1.M(fi).A(c0);
		System.out.println("ftest2 = "+ftest2);
		Function ftest3 = c1.D(fi).A(c3);
		System.out.println("ftest3 = "+ftest3);
		
		//f(x) = c1*l1(x)+c2
		FLinear1D l1 = new FLinear1D(1,5,2,20);
		System.out.println(l1);

		Function l2 = c1.M(l1).A(c2);
		System.out.println(l2);
		Variable var = new Variable(0.5);
		System.out.println(l2.value(var));
		
		Function ff = c1.D(fi).A(c3);
		Function ff_dx = ff._d("x");
		Function ff_dxx = ff_dx._d("x");
		System.out.println("ff = "+ff);
		System.out.println("ff_dx = "+ff_dx);
		System.out.println("ff_dxx = "+ff_dxx);
		System.out.println(ff_dx.value(var));
		
		//f(x) = 6*x^3 + 5*x^2 + 4*x + 3 多项式的导数
		List<Double> coef = new ArrayList<Double>();
		coef.add(3.0);
		coef.add(4.0);
		coef.add(5.0);
		coef.add(6.0);
		FPolynomial1D fp = new FPolynomial1D(coef);
		System.out.println(fp.value(new Variable(2.0)));
		Function fp_x2 = fp._d("x")._d("x");
		System.out.println(fp_x2.value(new Variable(2.0)));
		Function fp_x2d = (Function)fp_x2;
		Function fp_x3 = fp_x2d._d("X");
		System.out.println(fp_x3.value(new Variable(3.0)));
		
		Function power = FMath.pow(c2, c3);
		System.out.println(power+"="+power.value(null));
		
	}
	
	public static void testOperation() {
		Function fx = FX.fx;
		//f(x)=2*x+3
		Function f1 = new FAxpb(2.0,0.0);
		Function f2 = new FAxpb(2.0,3.0);
		Function f3 = new FAxpb(0.0,3.0);
		Function f4 = new FAxpb(0.0,0.0);
		System.out.println("f(x)="+f2);
		System.out.println(fx.M(f1));
		System.out.println(fx.M(f2));
		System.out.println(fx.M(f3));
		System.out.println(fx.M(f4));
		
		
	}
	
	public static void severalVariableFunctions() {
		Function fx = FX.fx;
		Function fy = FX.fy;
		Function f = FC.c(0.25).M(FC.C1.S(fx)).M(FC.C1.S(fy));
		System.out.println(f);
		Variable v = new Variable("x",0.5).set("y", 0.5);
		System.out.println(f.value(v));
	}
	
	
	public static void testDuDx() {
	    MeshReader reader = new MeshReader("rectangle.grd");
	    Mesh mesh = reader.read2DMesh();
	    mesh.computeNodeBelongsToElements();
	    
	    Vector v = Tools.function2vector(mesh, FX.fx.M(FX.fy));
	    
	    Tools.plotVector(mesh, "testCase", "v.dat", v);
	    DuDx vx = new DuDx(mesh,new Vector2Function(v, mesh, "x","y"),"y");
	    Tools.plotFunction(mesh, "testCase", "v_x.dat", vx);
	    
	    Vector[] a = new Vector[4];
	    for(int i=0;i<4;i++)
	    	a[i] = new SparseVectorHashMap(mesh.getNodeList().size());
	    ElementList eList = mesh.getElementList();
	    for(int i=1;i<=eList.size();i++) {
	    	Element e  = eList.at(i);
	    	NodeList nList = e.nodes;
	    	System.out.print(e);
	    	for(int j=1;j<=nList.size();j++) {
	    		Node node = nList.at(j);
		    	Variable vv = new Variable();
		    	vv.setIndex(node.globalIndex);
		    	vv.setElement(e);
		    	double dvx = vx.value(vv);
		    	System.out.print(" "+dvx+" ");
	    	}
	    	System.out.println();
	    }
	    
	    //D(x*y)/Dy=x 
	    //On element 1: [-3,-2.3]*[2.3,3]
	    double d = vx.value(new Variable("x",-2.8).set("y",2.9).setElement(eList.at(1)));
	    System.out.println(d);//-2.8
	    
	    
	    
	    
	    
	}
	
	public static void testFunctionExpression() {
		List<String> varNames = new ArrayList<String>();
		varNames.add("x");
		varNames.add("y");
		
		VectorFunction f = new AbstractVectorFunction(3,"x","y") {
			protected Function[] data = new Function[dim];
			@Override
			public void set(int index, Function value) {
				data[index] = value;
			}
			@Override
			public Function get(int index) {
				return data[index];
			}
		};
		System.out.println(f);
		System.out.println(f.getExpression());
		f.setFName("phi(x,y)");
		System.out.println(f);
		System.out.println(f.getExpression());
	}
	
	public static void testLinearCombination() {
		double[] cs = new double[3];
		Function[] fs = new Function[3];
		cs[0] = 1;
		cs[1] = 2;
		cs[2] = 3;
		fs[0] = FX.fx;
		fs[1] = FX.fy;
		fs[2] = FX.fz;
		Function lc = FMath.linearCombination(cs, fs);
		System.out.println(lc);
		Function lcd = lc._d("x");
		System.out.println(lcd);
		VariableArray va = new VariableArray();
		double[] x = {10, 10};
		double[] y = {20, 20};
		double[] z = {30, 30};
		va.set("x", x);
		va.set("y", y);
		va.set("z", z);
		
		printArray(lc.valueArray(va, null));
		printArray(lcd.valueArray(va, null));
		
	}
	
	public static void printArray(double[] ary) {
		for(int i=0;i<ary.length;i++)
			System.out.print(ary[i]+" ");
		System.out.println();
	}
	
	
	public static void testForScala() {
	    //f1 = x^2 + 2x + 2
	    Function f1 = X.M(X).A( X.A(1.0).M(2.0) );
	    System.out.println(f1);
	    System.out.println(f1.value(new Variable(2.0)));
	    System.out.println(f1._d("x"));
	    System.out.println(f1._d("x").value(new Variable(3.0)));

	    //f2 = 3x^2 + 4y + 1
	    Function f2 = new AbstractFunction("x","y") {
	    	@Override
	    	public double value(Variable v) {
		          double x = v.get("x");
		          double y = v.get("y");
		          return 3*x*x + 4*y + 1; 
	    	}
	    	@Override
	    	public Function _d(String vn) {
	    		if (vn == "x") return C(6).M(X);
	    		else if(vn == "y") return C(4);
	    		else throw new FutureyeException("variable name="+vn);
	    	}
	    }.setFName("f(x,y) = 3x^2+4y+1");
	    
	    System.out.println(f2);
	    Variable vv = new Variable("x",2.0).set("y",3.0);
	    System.out.println(f2.value(vv));
	    
	    System.out.println(f2._d("x"));
	    System.out.println(f2._d("y"));
	    //System.out.println(f2._d("z")); //Exception
	    
	    //f3 = (x+y, x*y, 1.0)
	    VectorFunction f3 = new SpaceVectorFunction( 
	    		new AbstractFunction("x","y") {
	    	    	@Override
	    	    	public double value(Variable v) {
	    		          return v.get("x") + v.get("y"); 
	    	    	}
	    	    }.setFName("x+y"),
	    	    X.M(Y),
	    	    C1
	        );
	    System.out.println(f3);
	    System.out.println(f3.value(new Variable("x",2.0).set("y",2.0)));
	}
	
	public static void main(String[] args) {
//		test();
//		constantTest();
//		testOperation();
//		
//		severalVariableFunctions();
		//testDuDx();
//		testFunctionExpression();
		
//		testLinearCombination();
		
	FAxpb f = new FAxpb(2.0,1.0);
	double v = f.value(new Variable("x",3.0));
	System.out.println(v);
    Function df = f._d("x");
    System.out.println(df);
	
		
		testForScala();
	}

		
}
