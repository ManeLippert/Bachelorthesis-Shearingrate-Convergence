# Einführung in Plasmaphysik

#### Vorlesung von Wolfgang Suttrop im Master Physik an der Universität Bayreuth

## 1. Einleitung

Plasmen, ionisiertes Gas, unterscheiden sich in ihren Eigenschaften wesentlich von einem neutralen Gas und werden daher auch der "vierte Aggregatszustand" genannt. Der überwiegende Anteil der sichtbaren Materie ist im Plasmazustand: Alle Sterne, interstellares Gas, u.v.m. Die langreichweitigen elektromagnetischen Kräfte zwischen den geladenen Plasmateilchen verursachen eine Vielzahl von kollektiven Effekten in einem Plasma.

In der ersten Vorlesung betrachten wir:
- Beispiele für Plasmen auf sehr unterschiedlichen Größenskalen
- Räumliche und zeitliche Skalen für die Quasineutralität eines Plasmas
- Bedingung für kollektives Verhalten
- Zustandsgrenzen für ideale, nicht-relativistische und nicht-entartete Plasmen

[Vorlesung](/literatur/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/01-Einleitung/01-Einleitung.pdf)

## 2. Gasentladung

Gasentladungen sind eine Sammelbezeichnung für verschiedene technische Plasmen die durch elektrischen Stromfluss erzeugt und aufrechterhalten werden. Wir betrachten die grundlegenden Mechanismen für die drei Haupttypen von Gasentladungen:
- Unselbständige Gasentladungen (z.B. Geiger-Müller Zähler)
- Glimmentladungen  (z.B. Niederdruck-Leuchtmittel)
- Bogenentladungen  (z.B. Hochdruck-Leuchtmittel, Plasma-Oberflächenbehandlung)

In der ersten Vorlesung betrachten wir:
- Grundlagen der Ladungsträger-Erzeugung und Vervielfachung
- Elektrischen Durchbruch
- Glimmentladungen 
- Korona-Entladungen 
- Bogen-Entladungen

[Vorlesung](/literatur/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/02-Gasentladungen.pdf)

### Simulation

[xpdp1](/literatur/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/xpdp1/)


## 3. Ionisationsgrad des Plasmas

Wie stark ist das Plasma ionisiert? Wir kennen schwach ionisierte Plasmen (z.B. Gasentladungen bei kleinen Strömen) und stark ionisierte Plasmen (z.B. Sterninneres oder Fusionsplasmen). Bei hoher Rate der Ionisations- und Rekombinationsprozesse kann man die Besetzung der verschiedenen Teilchenzustände oft wie in einem thermischen Gleichgewicht gemäß den Gesetzen der statistischen Physik durch eine gemeinsame Temperatur beschreiben, obwohl ein Plasma eigentlich fast immer kein abgeschlossenes thermodynamisches System darstellt. Bei kleinen Wechselwirkungsraten genügt es oft, die jeweils stärksten Ionisations- und Rekombinationsprozesse zu bilanzieren (Ratengleichgewicht).

[Vorlesung](/literatur/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/03-Ionisationsgrad/03-Ionisationsgrad.pdf)

## 04 Coulomb-Stöße

In vielen Plasmen dominiert die Wechselwirkung durch die langreichweitige Coulomb-Wechselwirkung über (kurzreichweitige) atomare Anregungs- oder Ionisationsstöße. Coulomb-Stöße sind elastisch, d.h. erhalten Gesamtenergie und Gesamtimpuls der Stoßpartner, da keine inneren atomaren Anregungen erfolgen. Die Abbremsung von Ladungsträgern, Ursache für den elektrischen Widerstand des Plasmas,  erfolgt durch Ablenkung der Bewegungsrichtung. Aufgrund der langen Reichweite des elektrostatischen Potenzials dominieren Stöße mit kleinen Ablenkwinkeln. Für ein reines Coulomb-Potenzial divergiert der totale Wirkungsquerschnitt (wie aus der klassischen Mechanik bekannt) und die Abbremskraft. Für ein abgeschirmtes Potenzial wie es im quasineutralen Plasma durchweg vorliegt, können wir in üblicher Näherung die Kraftwirkung außerhalb einer Abschirmlänge (Debye-Länge) vernachlässigen und so eine effektive Stossfrequenz und damit den elektrischen Widerstand des Plasmas berechnen.

[Vorlesung](/literatur/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/04-Coulombstoesse/04-Coulombstoesse.pdf)

### Python

[coulomb.py](/literatur/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/04-Coulombstoesse/coulomb.py), 
[coulomb_animiert](/literatur/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/04-Coulombstoesse/coulomb_anim1.py),
[coulomb_statisch](/literatur/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/04-Coulombstoesse/coulomb_stat.py)