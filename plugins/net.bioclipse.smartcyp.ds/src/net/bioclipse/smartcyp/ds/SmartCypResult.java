package net.bioclipse.smartcyp.ds;

import java.util.Map;

import net.bioclipse.ds.model.result.AtomResultMatch;

import org.openscience.cdk.renderer.generators.IGeneratorParameter;

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
