package net.bioclipse.smartcyp.ds;

import java.util.Map;

import org.openscience.cdk.renderer.generators.IGeneratorParameter;

import net.bioclipse.cdk.domain.ICDKMolecule;
import net.bioclipse.ds.model.result.AtomResultMatch;
import net.bioclipse.ds.model.result.BlueRedColorScaleGenerator;
import net.bioclipse.ds.model.result.SimpleResult;

/**
 * 
 * @author ola
 *
 */
public class SmartCypResult extends AtomResultMatch {

	public SmartCypResult(String name, int classification) {
		super(name, classification);
	}
	
	@Override
	public Class<? extends IGeneratorParameter<Boolean>> getGeneratorVisibility() {
    	return (Class<? extends IGeneratorParameter<Boolean>>)SmartcypGenerator.Visibility.class;
    }

    @Override
    public Class<? extends IGeneratorParameter<Map<Integer, Number>>> getGeneratorAtomMap() {
    	return (Class<? extends IGeneratorParameter<Map<Integer, Number>>>)SmartcypGenerator.AtomMap.class;
    }
}
