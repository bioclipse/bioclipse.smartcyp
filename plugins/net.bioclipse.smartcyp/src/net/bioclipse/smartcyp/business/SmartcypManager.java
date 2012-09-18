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

import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import net.bioclipse.cdk.business.Activator;
import net.bioclipse.cdk.business.ICDKManager;
import net.bioclipse.cdk.domain.CDKMolecule;
import net.bioclipse.cdk.domain.ICDKMolecule;
import net.bioclipse.core.business.BioclipseException;
import net.bioclipse.core.domain.IMolecule;
import net.bioclipse.managers.business.IBioclipseManager;

import org.apache.log4j.Logger;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.AtomNumberGenerator;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.IGenerator;
import org.openscience.cdk.renderer.generators.RingGenerator;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;

import smartcyp.GenerateImages;
import smartcyp.MoleculeKU;
import smartcyp.SMARTSnEnergiesTable;
import smartcyp.rankedlabelgenerator;
import smartcyp.rankedlabelgenerator2C9;
import smartcyp.rankedlabelgenerator2D6;
import smartcyp.MoleculeKU.SMARTCYP_PROPERTY;

public class SmartcypManager implements IBioclipseManager {

    private static final Logger logger = Logger.getLogger(SmartcypManager.class);

	int WIDTH = 800;
	int HEIGHT = 400;
	Rectangle drawArea = new Rectangle(WIDTH, HEIGHT);
	Image image = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);

    
    /**
     * Gives a short one word name of the manager used as variable name when
     * scripting.
     */
    public String getManagerName() {
        return "smartcyp";
    }
    
    public ICDKMolecule predictSOM(IMolecule mol) throws BioclipseException, CloneNotSupportedException, CDKException {
    	
    	ICDKManager cdk = Activator.getDefault().getJavaCDKManager();
    	
    	ICDKMolecule cdkmol = cdk.asCDKMolecule(mol);
    	IAtomContainer ac = cdkmol.getAtomContainer();
    	SMARTSnEnergiesTable smartsEnergies = new SMARTSnEnergiesTable();
    	
    	MoleculeKU moleculeKU = new MoleculeKU(ac, smartsEnergies.getSMARTSnEnergiesTable());	
		moleculeKU.setProperty(CDKConstants.TITLE, ac.getProperty(CDKConstants.TITLE));
		moleculeKU.setProperties(ac.getProperties());

		moleculeKU.assignAtomEnergies(smartsEnergies.getSMARTSnEnergiesTable());	

		//System.out.println("\n ************** Calculating shortest distance to protonated amine **************");
		moleculeKU.calculateDist2ProtAmine();

		//System.out.println("\n ************** Calculating shortest distance to carboxylic acid **************");
		moleculeKU.calculateDist2CarboxylicAcid();
		
		//System.out.println("\n ************** Calculating Span2End**************");
		moleculeKU.calculateSpan2End();
		
		//System.out.println("\n ************** Calculating Accessabilities and Atom Scores**************");
		moleculeKU.calculateAtomAccessabilities();
		moleculeKU.calculateAtomScores();
		moleculeKU.calculate2D6AtomScores();
		moleculeKU.calculate2C9AtomScores();

		
		//System.out.println("\n ************** Identifying, sorting and ranking C, N, P and S atoms **************");
		moleculeKU.sortAtoms();
		moleculeKU.rankAtoms();
		moleculeKU.sortAtoms2D6();
		moleculeKU.rankAtoms2D6();
		moleculeKU.sortAtoms2C9();
		moleculeKU.rankAtoms2C9();
		
		//Debug out properties
//		System.out.println("SmartCyp properties per atom: ");
//		for(IAtom atom : moleculeKU.atoms() ){
//			System.out.println("Atom: " + moleculeKU.getAtomNumber(atom) +" has Ranking: " + SMARTCYP_PROPERTY.Ranking.get(atom));
//			for (Object prop : atom.getProperties().keySet()){
//				String key = (String) prop;
//				Object val = atom.getProperties().get(prop);
//				System.out.println("    " + key + " -- " + val);
//			}
//		}

//		for(IAtom atom : moleculeKU.getAtomsSortedByEnA() ){
//			System.out.println("Sorted atom: " + moleculeKU.getAtomNumber(atom));
//				Object val = atom.getProperty("Score");
//				System.out.println("    Score:" + val);
//		}

		return new CDKMolecule(moleculeKU);
		
/*
		System.out.println("Atom,Ranking,Score,Energy,Relative Span,2D6ranking,2D6score,Span2End,N+Dist,2C9ranking,2C9score,COODist");
		DecimalFormat twoDecimalFormat = new DecimalFormat("#.##");

		// Iterate Atoms
		for(int atomIndex = 0; atomIndex < moleculeKU.getAtomCount()  ; atomIndex++ ){
			
			IAtom currentAtom = (Atom) moleculeKU.getAtom(atomIndex);

			// Match atom symbol
			String currentAtomType = currentAtom.getSymbol();
			if(currentAtomType.equals("C") || currentAtomType.equals("N") || currentAtomType.equals("P") || currentAtomType.equals("S")) {
				
				System.out.print(currentAtom.getSymbol() + "."+ currentAtom.getID() + "," + SMARTCYP_PROPERTY.Ranking.get(currentAtom) + ",");				
				if(SMARTCYP_PROPERTY.Score.get(currentAtom) != null) 
					System.out.print(twoDecimalFormat.format(SMARTCYP_PROPERTY.Score.get(currentAtom)) + "," + SMARTCYP_PROPERTY.Energy.get(currentAtom));
				else System.out.print("999,999");
				System.out.print("," + twoDecimalFormat.format(SMARTCYP_PROPERTY.Accessibility.get(currentAtom)));
				if(SMARTCYP_PROPERTY.Score2D6.get(currentAtom) != null) {
					System.out.print("," + SMARTCYP_PROPERTY.Ranking2D6.get(currentAtom));
					System.out.print("," + twoDecimalFormat.format(SMARTCYP_PROPERTY.Score2D6.get(currentAtom)));
				}
				else System.out.print("999,999");
				System.out.print("," + twoDecimalFormat.format(SMARTCYP_PROPERTY.Span2End.get(currentAtom)));
				if(SMARTCYP_PROPERTY.Dist2ProtAmine.get(currentAtom) != null)
					System.out.print("," + twoDecimalFormat.format(SMARTCYP_PROPERTY.Dist2ProtAmine.get(currentAtom)));
				else System.out.print(",0");
				if(SMARTCYP_PROPERTY.Score2C9.get(currentAtom) != null) {
					System.out.print("," + SMARTCYP_PROPERTY.Ranking2C9.get(currentAtom));
					System.out.print("," + twoDecimalFormat.format(SMARTCYP_PROPERTY.Score2C9.get(currentAtom)));
				}
				else System.out.print("999,999");
				if(SMARTCYP_PROPERTY.Dist2CarboxylicAcid.get(currentAtom) != null)
					System.out.print("," + twoDecimalFormat.format(SMARTCYP_PROPERTY.Dist2CarboxylicAcid.get(currentAtom)));
				else System.out.print(",0");
				System.out.print("\n");
			}
		}
		
		
		// Generate 2D coordinates for moleculeKU
		ac = GenerateImages.generate2Dcoordinates(ac);

		// Generators make the image elements
		List<IGenerator<IAtomContainer>> generators = new ArrayList<IGenerator<IAtomContainer>>();
		generators.add(new BasicSceneGenerator());
		generators.add(new RingGenerator());
		generators.add(new BasicAtomGenerator());
		generators.add(new AtomNumberGenerator());
		
		String OutputDir = "/tmp/";

		// The renderer renders the picture
		//AWTFontManager fontman = new AWTFontManager();
		//renderer = new Renderer(generators, fontman, false);
		AtomContainerRenderer renderer =
			  new AtomContainerRenderer(generators, new AWTFontManager());
		renderer.setup(ac, drawArea);
		
		// Set layout of molecule
		// This method is not used because default layout looks ok
		// this.setMoleculeLayout(r2dm);


		// Write 2 types of images with: 1) heteroatoms and 2) atom Numbers
		renderer.getRenderer2DModel().set(AtomNumberGenerator.ColorByType.class, true);
		this.paintAndWriteMolecule(renderer, ac, "atomNumbers", OutputDir);
		generators.removeAll(generators);
		generators.add(new BasicSceneGenerator());
		generators.add(new rankedlabelgenerator());
		generators.add(new RingGenerator());
		generators.add(new BasicAtomGenerator());
		this.paintAndWriteMolecule(renderer, ac, "heteroAtoms", OutputDir);	
		generators.removeAll(generators);
		generators.add(new BasicSceneGenerator());
		generators.add(new rankedlabelgenerator2D6());
		generators.add(new RingGenerator());
		generators.add(new BasicAtomGenerator());
		this.paintAndWriteMolecule(renderer, ac, "heteroAtoms2D6", OutputDir);
		generators.removeAll(generators);
		generators.add(new BasicSceneGenerator());
		generators.add(new rankedlabelgenerator2C9());
		generators.add(new RingGenerator());
		generators.add(new BasicAtomGenerator());
		this.paintAndWriteMolecule(renderer, ac, "heteroAtoms2C9", OutputDir);
*/
	}
    
	public void paintAndWriteMolecule(AtomContainerRenderer renderer, IAtomContainer iAtomContainer, String nameBase, String outputdir){


		// Paint background
		Graphics2D g2 = (Graphics2D)image.getGraphics();
		//	g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, WIDTH, HEIGHT);

		// the paint method also needs a toolkit-specific renderer
		renderer.paint(iAtomContainer, new AWTDrawVisitor(g2), drawArea, true);
		
		String moleculeID = iAtomContainer.getID();
		//System.out.println(moleculeID);
		String fileName = outputdir + "smartcyp_images" + File.separator + "molecule" + "_" + nameBase + ".png";
		
		System.out.println(fileName);
		try {ImageIO.write((RenderedImage)image, "PNG", new File(fileName));} 
		catch (IOException e) {
			e.printStackTrace();
			System.out.println("Molecule images could not be written to file. " +
			"If you have not already you need to create the directory 'smartcyp_images' in which the images are to be written");
		}

	}
    
}
