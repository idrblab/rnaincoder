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

@XmlRootElement(name = "loop")
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(propOrder = { "name", "strand" })
@NoArgsConstructor
@EqualsAndHashCode(exclude = { "strand" })
public class Loop implements DpParameter {
    @XmlElement
    @Getter
    @Setter
    private String name;

    @XmlElement
    @Getter
    @Setter
    private DpStrand strand;

    public Loop(final String[] elements) {
        this.name = elements[0];
        this.strand = new DpStrand(elements[1], elements[2]);
    }

    @Override
    public String toString() {
        return new StringBuilder("(\"").append(name).append("\",").append(strand.toString()).append(')')
                .toString();
    }
}
