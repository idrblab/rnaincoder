package edu.put.ma.templates;

public interface ResidueTemplate {

    String getOneResidueName();

    String getThreeResidueName();

    void addAtomTemplate(String atomNames, String templateAtomName);

    boolean isEqualResidueName(String presName);

    boolean isProperAtomName(String atomName);

    String getTemplateAtomName(String atomName);

    String getOriginalAtomName(String atomName);

    boolean isAppropriateResName(String resName);
}
