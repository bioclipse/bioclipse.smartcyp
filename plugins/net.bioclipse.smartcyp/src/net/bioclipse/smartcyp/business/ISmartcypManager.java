/*******************************************************************************
 * Copyright (c) 2012  Ola Spjuth <ola.spjuth@gmail.com>
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contact: http://www.bioclipse.net/
 ******************************************************************************/
package net.bioclipse.smartcyp.business;

import org.openscience.cdk.exception.CDKException;

import net.bioclipse.core.PublishedClass;
import net.bioclipse.core.PublishedMethod;
import net.bioclipse.core.Recorded;
import net.bioclipse.core.business.BioclipseException;
import net.bioclipse.core.domain.IMolecule;
import net.bioclipse.managers.business.IBioclipseManager;

@PublishedClass(
    value="TODO: Describe the manager here."
)
public interface ISmartcypManager extends IBioclipseManager {


    @Recorded
    @PublishedMethod(
    	params="IMolecule mol",
        methodSummary=
            "Predict Site-Of-Metabolism for a molecule using SmartCyp."
    )
    public void predictSOM(IMolecule mol) throws BioclipseException, CloneNotSupportedException, CDKException;
	
}
