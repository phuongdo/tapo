package org.cnrs.crbm.lib.services;

import org.eclipse.jetty.util.log.Logger;
import org.slf4j.LoggerFactory;

public class ServerLogger implements Logger {

	private org.slf4j.Logger _logger;

	public ServerLogger() {
		this("org.eclipse.jetty.util.log");
	}

	public ServerLogger(String name) {
		this._logger = LoggerFactory.getLogger(name);
	}

	public String getName() {
		return "Server Logging";
	}

	public void warn(String msg, Object... args) {
		this._logger.warn(msg, args);
	}

	public void warn(Throwable thrown) {
	}

	public void warn(String msg, Throwable thrown) {
	}

	public void info(String msg, Object... args) {
		this._logger.info(msg, args);
	}

	public void info(Throwable thrown) {
		info("", thrown);
	}

	public void info(String msg, Throwable thrown) {
		this._logger.info(msg, thrown);
	}

	public boolean isDebugEnabled() {
		return false;
	}

	public void setDebugEnabled(boolean enabled) {
	}

	public void debug(String msg, Object... args) {
	}

	public void debug(Throwable thrown) {
	}

	public void debug(String msg, Throwable thrown) {
	}

	public Logger getLogger(String name) {
		return this;
	}

	public void ignore(Throwable ignored) {
	}
}