//package org.cnrs.crbm.lib.trsfinder;
//
//
//import org.biojava.nbio.structure.StructureException;
//import org.cnrs.crbm.lib.repeats.module.Pro2Vector;
//
//@Deprecated
//public class VectorTReksFinder extends Finder {
//
//	public VectorTReksFinder(Features features) {
//		super(features);
//
//	}
//
//	@Override
//	public void findRepeat(Features features) {
//		try {
//			Pro2Vector pro2Vector = new Pro2Vector();
//			this.repeats = pro2Vector.getRepeatsTReks(
//					this.features.getVectors(), this.features.getAtoms());
//
//		} catch (StructureException e) {
//
//			e.printStackTrace();
//		}
//	}
//
//}
