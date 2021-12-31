package org.cnrs.crbm.lib.conf;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Properties;

/**
 * 
 * @author Phuong Read configure from file.
 */
public class ConfigUtil {
	private static final String CONFIG_FILE = Enviroment.getHome()
			+ "conf/struct-lib.cfg";
	private static ConfigUtil config = null;

	Reader reader = null;
	Properties prop = new Properties();

	/**
	 * Return configuration
	 */
	protected ConfigUtil() {

		try {
			reader = new FileReader(CONFIG_FILE);
			prop.load(reader);

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * 
	 * @return appconfig �ntance
	 */
	public static ConfigUtil getInstance() {

		if (config == null) {
			config = new ConfigUtil();
		}

		return config;
	}

	public String getProperty(String name) {
		return prop.getProperty(name);
	}

}
