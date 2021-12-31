package org.cnrs.crbm.lib.sadb;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ConformationMatcher {

    // example comes from 3iox_A
    public static final String EXAMPLE_TEST = "xbpaafaaaaaafaaaaaaaaaaaaafaaaaaafaaaaaglabbeaaaaafaaaaaagaaafaaaaaafaaaaaafaaaglxx";

	public static void main(String[] args) {

        String regex = ".*[agfh]{40,}.*";
        // System.out.println(EXAMPLE_TEST.matches(".*" + regex + ".*"));

        //System.out.println(isLongHelix(EXAMPLE_TEST));

		Pattern pattern = Pattern.compile(regex);
		// in case you would like to ignore case sensitivity,
		// you could use this statement:
		// Pattern pattern = Pattern.compile("\\s+", Pattern.CASE_INSENSITIVE);
		Matcher matcher = pattern.matcher(EXAMPLE_TEST);

		// check all occurance
		while (matcher.find()) {

            System.out.println(matcher.start() + "-" + matcher.end());
            //System.out.println(matcher.group());
        }

	}

	/**
	 * with h patterns >40
	 * 
	 * @param s
	 * @return
	 */
	public static boolean isLongHelix(String s) {
		String regex = ".*[agfh]{40,}.*";
		return s.matches(regex);
	}

	/**
	 * - b/p > 20 ( pdb: 1CAG ) => class II ( collagen) - P > 20 => class II
	 * 
	 * 
	 * @param s
	 * @return
	 */
	public static boolean isClass2(String s) {
		String regex = ".*[bpvs]{20,}.*";
		return s.matches(regex);
	}

	public static boolean isClass2OnlyP(String s) {
		String regex = ".*[ps]{15,}.*";
		return s.matches(regex);
	}

	public static boolean isClass1(String s) {
		String regex = ".*[bv]{15,}.*";
		return s.matches(regex);
	}

	/**
	 * with HBHBHB >40
	 * 
	 * @param s
	 * @return
	 */
	public static boolean isRossmanfoldPattern(String s) {
		String regex1 = ".*(HB){3,}.*";
		String regex2 = ".*(HB){3,}.*";
		return s.matches(regex1) || s.matches(regex2);
	}
}
