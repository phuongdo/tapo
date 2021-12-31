package org.cnrs.crbm.lib.utils;

import java.text.*;
import java.util.Locale;

public class NumberFormatUtils {

	private static final DecimalFormat df;
	static DecimalFormatSymbols decimalFormatSymbols = new DecimalFormatSymbols();

	static {
		decimalFormatSymbols.setDecimalSeparator('.');
		decimalFormatSymbols.setGroupingSeparator(',');

		df = new DecimalFormat("###0.000", decimalFormatSymbols);
	}// end static initialiser

	public static String format(double number) {
		return df.format(number);
	}

}// end class
