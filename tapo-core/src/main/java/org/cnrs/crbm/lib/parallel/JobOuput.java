package org.cnrs.crbm.lib.parallel;

import java.util.HashMap;
import java.util.Map;

public class JobOuput {

	Map<String, String> proteins = new HashMap<String, String>();

	public Map<String, String> getProteins() {
		return proteins;
	}

	public void setProteins(Map<String, String> proteins) {
		this.proteins = proteins;
	}

	public void join(JobOuput jobOuput) {
		proteins.putAll(jobOuput.getProteins());
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		for (Map.Entry<String, String> entry : this.proteins.entrySet()) {
			builder.append(entry.getValue() + "\n");
		}
		return builder.toString();
	}
}
