package edu.put.ma.model.services;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@XmlRootElement(name = "helix")
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(propOrder = { "name", "strand5", "strand3" })
@NoArgsConstructor
@EqualsAndHashCode(exclude = { "strand5", "strand3" })
public class Helix implements DpParameter {
    @XmlElement
    @Getter
    @Setter
    private String name;

    @XmlElement
    @Getter
    @Setter
    private DpStrand strand5;

    @XmlElement
    @Getter
    @Setter
    private DpStrand strand3;

    public Helix(final String[] elements) {
        this.name = elements[0];
        this.strand5 = new DpStrand(elements[1], elements[2]);
        this.strand3 = new DpStrand(elements[3], elements[4]);
    }

    @Override
    public String toString() {
        return new StringBuilder("(\"").append(name).append("\",").append(strand5).append(',')
                .append(strand3).append(')').toString();
    }
}
