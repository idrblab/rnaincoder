package edu.put.ma.templates;

public interface Templates {

    String getTemplatesString();

    void load(String dictionaryString);

    String getThreeLetterResidueNameByStructuralName(String resName);

    String getOneLetterResidueNameByStructuralName(String resName);

    String getTemplateAtomName(String resName, String atomName);

    String getOriginalAtomName(String resName, String atomName);

    boolean isProperAtomName(String residueName, String atomName);

    boolean isInitiated();
}
