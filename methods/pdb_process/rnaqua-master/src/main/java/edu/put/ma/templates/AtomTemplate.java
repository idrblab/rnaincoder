package edu.put.ma.templates;

public interface AtomTemplate {

    String getTemplateAtomName();
    
    String getOriginalAtomName();
    
    boolean isProperAtomName(String patomName);
}
