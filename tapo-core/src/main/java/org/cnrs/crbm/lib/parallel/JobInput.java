package org.cnrs.crbm.lib.parallel;

import java.util.ArrayList;
import java.util.List;

import org.cnrs.crbm.lib.io.Row;

public class JobInput {

	List<Row> proteins = new ArrayList<Row>();

	public List<Row> getProteins() {
		return proteins;
	}

	public void setProteins(List<Row> proteins) {
		this.proteins = proteins;
	}

}
