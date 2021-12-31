package org.cnrs.crbm.lib.domain;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.domain.pdp.CutDomain;
import org.biojava.nbio.structure.domain.pdp.Domain;
import org.biojava.nbio.structure.domain.pdp.Segment;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.List;

public class DemoDomainsplit {

	public static void main(String[] args) {

		DemoDomainsplit split = new DemoDomainsplit();

		// String pdbId = "3gly";
		String pdbId = "3iox";

		split.basicLoad(pdbId);

	}

	public void basicLoad(String pdbId) {

		try {

			// end of optional part

			Structure struc = PdbTools.getStructureFromLocalPdb(pdbId);

			// System.out.println("structure loaded: " + struc);

			CutDomain.verbose = false;
			Atom[] ca = StructureTools.getAtomCAArray(struc);
			List<Domain> domains = LocalProteinDomainParser.suggestDomains(ca);

			System.out.println("RESULTS: =====");
			for (Domain dom : domains) {
				System.out.println("DOMAIN:" + dom.getSize() + " "
						+ dom.getScore());
				List<Segment> segments = dom.getSegments();
				for (Segment s : segments) {
					System.out.println("   Segment: " + s);
					System.out
							.println(">> "
									+ ca[s.getFrom()].getGroup()
											.getResidueNumber()
									+ "-"
									+ ca[s.getTo()].getGroup()
											.getResidueNumber());

				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}