package edu.put.ma.model.input;

import org.apache.commons.cli.Options;

import edu.put.ma.Command;
import edu.put.ma.model.services.StructuresSet;

public interface CommonInputModel {

    boolean isInputInitializedProperly();

    Options constructSpecificOptions();

    void printHelp(String artifactId);

    String getInputModelString();

    String[] getArgs();

    StructuresSet getStructuresSet();

    String getOutputFilePath();

    String getServiceCode();

    Command getCommand();
    
    String getCommandName();

    StructuresSet constructStructuresSet();
}
