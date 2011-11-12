package edu.uta.futureye.core.intf;


import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ShapeFunction;

public interface WeakForm {
	static enum ItemType {Domain, Border};
	
	//--- Common approach to provide weak form to assembler-----
	//TODO 考虑直接传入DOF2011-11-6 (for upwind)
	void setShapeFunction(ShapeFunction trial, int trialDofLocalIndex,
			ShapeFunction test, int testDofLocalIndex);
	
	Function leftHandSide(Element e, ItemType itemType);
	Function rightHandSide(Element e, ItemType itemType);
	//----------------------------------------------------------
	
	//--- Fast approach to provide weak form to assembler-----
	/**
	 * Assemble element <code>e</code> here, instead of providing left hand side
	 * and right hand side.
	 * 
	 * @param e
	 * @param globalStiff (I/O): Global stiff matrix 
	 * @param globalLoad (I/O): Global load vector
	 *   
	 */
	void assembleElement(Element e, 
			Matrix globalStiff, Vector globalLoad);
	//--------------------------------------------------------
	
	/**
	 * Integrate on element <code>e</code>
	 * 
	 * @param e
	 * @param fun: LHS or RHS
	 * @return
	 */
	double integrate(Element e, Function fun);
	
}
