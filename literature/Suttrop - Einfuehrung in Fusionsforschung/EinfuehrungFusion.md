##### [Startpage](/README.md) &nbsp; &nbsp; &nbsp; &nbsp; [Journal](/journal/JOURNAL.md) &nbsp; &nbsp; &nbsp; &nbsp; Literature &nbsp; [[1]](/literature/Peeters%2C%20Rath%2C%20Buchholz%20-%20Gradient-driven%20flux-tube%20simulations%20of%20ion%20temperature%20gradient%20turbulence%20close%20to%20the%20non-linear%20threshold%20(Paper%2C%202016).pdf) [[2]](/literature/Peeters%2C%20Rath%2C%20Buchholz%20-%20Comparison%20of%20gradient%20and%20flux%20driven%20gyro-%0Akinetic%20turbulent%20transport%20(Paper%2C%202016).pdf) [[3]](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/EinfuehrungPlasma.md) [[4]](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/EinfuehrungFusion.md)
# Einführung in die Fusionsforschung

#### Vorlesung von [Dr. Wolfgang Suttrop](https://www.ipp.mpg.de/4123258/suttrop) im Master Physik an der Universität Bayreuth

## Inhaltsangabe
1.  [Energie aus Fusion](#1-energie-aus-fusion)
2.  [Plasmaphysik](#2-plasmaphysik)
3.  [Einschluss von Plasma im Magnetfeld](#3-einschluss-von-plasmen-im-magnetfeld)
4.  [Toroidale Kraftgleichgewichte](#4-toroidale-kraftgleichgewichte)
5.  [Driftbahnen im Torus](#5-driftbahnen-im-torus)
6.  [Neoklassischer Transport](#6-neoklassischer-transport)
7.  [Neoklassische Effekte](#7-neoklassische-effekte)
8.  [Plasma Diagnostik](#8-plasma-diagnostik)
9.  [Plasmarandschicht](#9-plasmarandschicht)
10. [Leistungsabfuhr](#10-leistungsabfuhr)
11. [Auslegung und Operationsgrenzen eines Fusionsreaktors](#11-auslegung-und-operationsgrenzen-eines-fusionsreaktors)
## 1. Energie aus Fusion

Einführung in das Thema der Vorlesung.
Wir rekapitulieren kurz die wichtigen Größen der Plasmaphysik.
Danach wenden wir uns den Voraussetzungen für die Energiegewinnung durch Kernfusion zu.
Das Lawson-Kriterium beschreibt, unter welchen Voraussetzungen sich ein Fusionsplasma selbst heizt (gegen Energieverluste) und netto Energie nach außen abgibt.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/01-Fusion/01-Fusion.pdf)


## 2. Plasmaphysik

Zusammenfassung der wichtigsten Resultate der Vorlesung "Einführung in die Plasmaphysik" mit Betonung der für die Fusionsforschung wichtigen Grundlagen. Dieses ist keine eigene Vorlesung, die Folien sollen als Kurzreferenz auf die benötigten Grundlagen der Vorlesung "Einführung in die Plasmaphysik" dienen.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/02-Plasma/02-Plasma.pdf)

[Komplette Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/EinfuehrungPlasma.md)


## 3. Einschluss von Plasmen im Magnetfeld

Plasmen -ionisierte Gase- bestehen aus geladenen Teilchen, die in elektrischen und magnetischen Feldern Kräfte erfahren. Wir betrachten Plasmen im sog. Flüssigkeitsbild, in dem Mittelwerte über die Geschwindigkeitsverteilung aller Plasmateilchen betrachtet werden. Unsere Plasmakonfiguration soll für lange Zeiten stationär vorliegen, daher müssen die treibenden Kräfte untereinander im Kraftgleichgewicht sein. Im einfachsten magnetischen Kraftgleichgewicht wird der im Plasma vorhandene Gradient des isotropen hydrostatischen Drucks durch die Lorentzkraft, die im Einflüssigkeits (MHD-) Modell die Form jxB annimmt, unterstützt. Wir betrachten grundlegende Eigenschaften von MHD-Kraftgleichgewichten und einfache, auch im Experiment zugängliche, Konfigurationen.

[Verlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/03-MagnetischerEinschluss/03-MagnetischerEinschluss.pdf)


## 4. Toroidale Kraftgleichgewichte

Einfache lineare Konfigurationen für den magnetischen Plasmaeinschluss haben den Nachteil, dass ein Experiment räumlich begrenzt ist und das Magnetfeld an den Enden die Wand schneiden muss. Da das Plasma entlang des Magnetfelds frei strömen kann, treten erhebliche Wärme- und Teilchen-Verluste auf. Die einfache Lösung besteht darin, das Magnetfeld toroidal in sich zu schliessen. Die Plasmawand wird dann nur durch (langsamen) Transport senkrecht zum Magnetfeld erreicht und das Plasma wird wesentlich besser eingeschlossen. Wir betrachten den einfachen Fall des axisymmetrischen Torus und leiten die sog. Grad-Shafranov-Schlüter-Gleichung ab, die den poloidalen magnetischen Fluss als Funktion von Radius und axialer Koordinate beschreibt (der toroidale Winkel ist wg. Axisymmetrie ignorabel). Als Beispiel für eine solche axisymmetrische toroidale Konfiguration betrachten wir den sog. Tokamak, der auch in der Fusionsforschung besonders erfolgreich ist. Anhand des Garchinger Tokamaks "ASDEX Upgrade" verdeutlichen wir uns den Aufbau des Experimentes und den Ablauf der darin untersuchten Plasmaentladungen.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/04-ToroidaleKonfigurationen/04-ToroidaleKonfigurationen.pdf)

### Simulation
[solovjov.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/04-ToroidaleKonfigurationen/solovjov.py)


## 5. Driftbahnen im Torus

Im toroidalen Magnetfeld gelingt der Plasma-Einschluss, wenn die vertikale Gyrozentrumsdrift (kombinierte Krümmungs- und <img src="https://render.githubusercontent.com/render/math?math={\color{white}\nabla\cdot\vec{B}}">-Drift) durch ein poloidales Magnetfeld kompensiert wird. Es bleiben Exkursionen der Gyrozentren weg von den magnetischen Flächen. Aufgrund der Inhomogenität des Magnetfelds (es ist innen im Torus höher als außen) ergibt sich sich ein magnetischer Spiegel, den nur Teilchen mit hoher Geschwindigkeit parallel <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{B}}"> umlaufen können und in dem Teilchen mit kleinerer Parallel-Geschwindigkeit auf der Niedrigfeldseite gefangen bleiben. Wir berechnen die Teilchenbahn beider Typen und finden dass gefangene Teilchen deutlich größere Exkursionen von der magnetischen Fläche unternehmen.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/05-Torusdrift/05-Torusdrift.pdf)

### Simulation
[tordrif_anim.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/05-Torusdrift/tordrift_anim.py)


## 6. Neoklassischer Transport

Wir betrachten grundsätzlich den Transport senkrecht zum Magnetfeld - im Teilchen- und im Flüssigkeitsbild, die beide im wesentlichen das gleiche Bild für die Teilchendiffusion darbieten. Im Torus (mit zwangsweise inhomogenem und gekrümmtem Magnetfeld) ergeben sich wichtige Modifikationen: Der "neoklassische Transport" ist  gegenüber dem Transport im homogenen Magnetfeld erheblich durch die Anwesenheit von Teilchendriften erhöht. Je nach Verhältnis von Stoßfrequenz und Transitfrequenz (Frequenz der poloidalen Gyrozentrums-Umläufe) unterscheidet man verschiedene Regimes: Im "Bananen"-Regime (Grenzfall kleiner Stossfrequenzen) wird der radiale Transport durch im poloidalen magnetischen Spiegel gefangene Teilchen dominiert, obwohl diese nur einen kleinen Bruchteil der gesamten Teilchen-Population ausmachen. Im "Pfirsch-Schlüter"-Regime (Grenzfall hoher Stossfrequenzen) wird der radiale Transport durch alle Teilchen getragen, ist aber wg. der Driften weiterhin gegenüber dem Fall homogenen Magnetfelds erhöht, wenn auch schwächer als im Bananen-Regime.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/06-NeoklassicherTransport/06-NeoklassischerTransport.pdf)


## 7. Neoklassische Effekte

Wir behandeln weitere neoklassische Effekte: Den "Ripple"-Transport durch nicht-axisymmetrisches Magnetfeld, den "Bootstrap"-Strom durch gefangene Teilchen bei Anwesenheit eines radialen Druckgradienten, und die <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{E}\times\vec{b}}"> sowie den "Ware-pinch" bei Anwesenheit eines toroidalen elektrischen Feldes (wie jenes das im Tokamak den toroidalen Plasmastrom antreibt).

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/07-NeoklassischeEffekte/07-NeoklassischeEffekte.pdf)


## 8. Plasma-Diagnostik

Zur Erforschung der Plasma-Eigenschaften sowie zur Steuerung und Regelung eines möglichen Fusionsreaktors ist die Messung vielen Kenngrößen erforderlich. Im heissen und dichten Plasma müssen Messungen i.w. "berührungslos", d.h. durch elektrische, magnetische oder optische (elektromagnetische) Effekte des Plasmas erfolgen. Wir diskutieren Messmethoden die für den magnetischen Einschluss wichtig sind: Magnetische Messungen, Messmethoden mit Millimeter-Wellen, Methoden der Plasmaspektroskopie (sichtbares und UV-Licht).

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/08-Diagnostik/08-Diagnostik.pdf)


## 9. Plasmarandschicht

Im brennenden Fusionsplasma wird die durch die erzeugten energetischen Alpha-Teilchen (He-Kerne) im Plasma abgegebene Wärme durch radialen Transport auf die das Plasma umgebende Wand abgegeben. Je nach Konfiguration können dabei sehr hohe Leistungsdichten auftreten, die das Wandmaterial gefährden können. Wir betrachten zunächst die Verhältnisse in der sog. elektrostatische Plasmarandschicht, eine nicht-neutrale dünne Schicht direkt an der Oberfläche der Wand, die Teilchen- und Wärmeflüsse auf die Wand reguliert. Danach untersuchen wir den Plasmarand, d.h. die Zone zwischen geschlossenen magnetischen Flächen und der elektrostatischen Randschicht, in der die Wärme- und Teilchenflüsse Gradienten der Temperatur- und Dichte bewirken. Je nach Plasmadichte unterscheidet man das sog. "Randschicht-begrenzte" und das sog. "Hoch-Recycling"-Regime. In Letzterem wird durch hohe Dichte direkt an der Wand die Temperatur klein gehalten, was besonders günstig für die Leistungsabfuhr ist.

[Vorlesung (1)](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/09-Plasmarandschicht/09-Randschicht1.pdf), [Vorlesung (2)](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/09-Plasmarandschicht/09-Randschicht2.pdf)

### Simulation

[randschicht.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/09-Plasmarandschicht/randschicht.py)

[xpdp1](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/xpdp1/) ([qmach2.inp](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/xpdp1/inp/qmach2.inp))

## 10. Leistungsabfuhr

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/10-Leistungsabfuhr/10-Leistungsabfuhr.pdf)


## 11. Auslegung und Operationsgrenzen eines Fusionsreaktors

Aufgrund einfacher Zusammenhänge können wir die wesentlichen Kenngrößen eines Fusionsreaktors bestimmen: Toroidales Magnetfeld, Plasmadichte, für das Brennen benötigte und tatsächlich erwartete Einschlusszeit. Letztere Bedingung ergibt eine Anforderung für die Mindestgröße eines energieliefernden Reaktors - wir erhalten ziemlich genau die Parameter des ersten "echten" Fusionsreaktors ITER, der gerade in Südfrankreich im Bau ist.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/11-Fusionsreaktor/11-Fusionsreaktor.pdf)

### Simulation
[lawson.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/11-Fusionsreaktor/lawson.py), [tokscaling.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Fusionsforschung/11-Fusionsreaktor/tokscaling.py)