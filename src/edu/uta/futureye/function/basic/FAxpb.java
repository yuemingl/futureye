package edu.uta.futureye.function.basic;

import java.util.List;

import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.util.Constant;

public class FAxpb extends AbstractFunction {
	protected double a;
	protected double b;
	//varName在构造的时候已经确定，以后不可以修改，但是varNames可以修改
	protected String varName;

	public FAxpb(double a, double b) {
		varNames.add(Constant.x);
		varName = Constant.x;
		this.a = a;
		this.b = b;
	}
	
	public FAxpb(String varName, double a, double b) {
		super(varName);
		this.varName = varName;
		this.a = a;
		this.b = b;
	}
	
	@Override
	public Function _d(String varName) {
		if(this.varName.equals(varName))
			return new FC(a);
		else
			return FC.c0;
	}

	@Override
	public double value(Variable v) {
		double rlt = a*v.get(varName)+b;
		return rlt;
	}
	
	@Override
	public int getOpOrder() {
		if(Double.compare(a, 0.0) == 0)
			return OP_ORDER0;
		if(Double.compare(b, 0.0) == 0)
			return OP_ORDER2;
		else
			return OP_ORDER3;
	}
	
	public String toString() {
		if(Double.compare(a, 1.0) == 0) {
			if(Double.compare(b, 0.0) == 0)
				return varName;
			else
				return varName+"+"+b;
		} else if(Double.compare(a, 0.0) == 0) {
				return b+"";
		} else if(Double.compare(b, 0.0) == 0) {
			return a+"*"+varName;
		}
		return a+"*"+varName+"+"+b;
	}
	
	/**
	 * varNames在由单个自变量表达式运算组合而成的多自变量函数情况下计算导数后可能被修改，
	 * 这种修改是允许的，例如：(x+1)*(y+1)关于x求导数后，原来(y+1)的varNames由[y]变为[x,y]
	 * 
	 */
	@Override
	public Function setVarNames(List<String> varNames) {
		this.varNames = varNames;
		return this;
	}
}
