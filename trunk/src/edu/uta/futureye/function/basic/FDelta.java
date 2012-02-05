package edu.uta.futureye.function.basic;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;

public class FDelta extends AbstractFunction {
	Variable x0;
	double eps;
	double amp;

	/**
	 * Delta(x) = amp * e^(-0.25*|x-x0|^2/eps ) / (2*sqrt(PI*eps))
	 * 
	 * @param x0
	 * @param eps
	 * @param amp
	 */
	public FDelta(Variable x0,double eps,double amp) {
		if(x0.getValues().size() == 1) {
			varNames.add("x");
		} else if(x0.getValues().size() == 2) {
			varNames.add("x");
			varNames.add("y");
		} else if(x0.getValues().size() == 3) {
			varNames.add("x");
			varNames.add("y");		
			varNames.add("z");		
		}
		this.x0 = x0;
		this.eps = eps;
		this.amp = amp;
	}
	
	@Override
	public double value(Variable x) {
		double d2 = 0.0;
		if(x0.getValues().size() == 1) {
			double dx = x.get("x")-x0.get("x");
			d2 = dx*dx;
		} else if(x0.getValues().size() == 2) {
			double dx = x.get("x")-x0.get("x");
			double dy = x.get("y")-x0.get("y");
			d2 = dx*dx+dy*dy;
		} else if(x0.getValues().size() == 3) {
			double dx = x.get("x")-x0.get("x");
			double dy = x.get("y")-x0.get("y");
			double dz = x.get("z")-x0.get("z");
			d2 = dx*dx+dy*dy+dz*dz;
		}
		return amp*Math.exp(-d2/eps/4.0) / (2*Math.sqrt(Math.PI*eps));
	}

	public double value() {
		return 0.0;
	}
	
	@Override
	public Function _d(String varName) {
		return null;
	}
	
	public String toString() {
		return "Delta";
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Variable x0 = new Variable();
		x0.set("x", 0.5);
		x0.set("y", 0.5);
		
		Variable x = new Variable();
		x.set("x", 0.6);
		x.set("y", 0.6);
		
		FDelta delta = new FDelta(x0,1e-4,1e5);
		System.out.println(delta.value(x));
	}

}
