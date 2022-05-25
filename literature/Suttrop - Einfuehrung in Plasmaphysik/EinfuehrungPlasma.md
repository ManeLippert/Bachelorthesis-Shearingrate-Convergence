##### [Startpage](/README.md) &nbsp; &nbsp; &nbsp; &nbsp; [Journal](/journal/JOURNAL.md) &nbsp; &nbsp; &nbsp; &nbsp; [Literature](/README.md#litarture) &nbsp; &nbsp; &nbsp; &nbsp; [GitHub](https://github.com/ManeLippert/Bachelorthesis-ZonalFlows)

# Einführung in Plasmaphysik

#### Vorlesung von [Dr. Wolfgang Suttrop](https://www.ipp.mpg.de/4123258/suttrop) im Master Physik an der Universität Bayreuth

## Inhaltsangabe
1.  [Einleitung](#1-einleitung)
2.  [Gasentladung](#2-gasentladung)
3.  [Ionisationsgrad des Plasmas](#3-ionisationsgrad-des-plasmas)
4.  [Coulomb-Stoesse](#4-coulomb-stoesse)
5.  [Einzelteilchen-Bewegung](#5-einzelteilchen-bewegung)
6.  [Hochfrequenz-Plasmaentladungen](#6-hochfrequenz-plasmaentladungen)
7.  [ Kinetische Beschreibung von Plasmen](#7-kinetische-beschreibung-von-plasmen)
8.  [Landau Dämpfung](#8-landau-dämpfung)
9.  [Flüssigkeitsbeschreibung von Plasma](#9-flüssigkeitsbeschreibung-von-plasma)
10. [Wellen im kalten Plasma](#10-wellen-im-kalten-plasma)
11. [Magnetohydrodynmaik](#11-magnetohydrodynamik)

## 1. Einleitung

Plasmen, ionisiertes Gas, unterscheiden sich in ihren Eigenschaften wesentlich von einem neutralen Gas und werden daher auch der "vierte Aggregatszustand" genannt. Der überwiegende Anteil der sichtbaren Materie ist im Plasmazustand: Alle Sterne, interstellares Gas, u.v.m. Die langreichweitigen elektromagnetischen Kräfte zwischen den geladenen Plasmateilchen verursachen eine Vielzahl von kollektiven Effekten in einem Plasma.

In der ersten Vorlesung betrachten wir:
- Beispiele für Plasmen auf sehr unterschiedlichen Größenskalen
- Räumliche und zeitliche Skalen für die Quasineutralität eines Plasmas
- Bedingung für kollektives Verhalten
- Zustandsgrenzen für ideale, nicht-relativistische und nicht-entartete Plasmen

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/01-Einleitung/01-Einleitung.pdf)

## 2. Gasentladung

Gasentladungen sind eine Sammelbezeichnung für verschiedene technische Plasmen die durch elektrischen Stromfluss erzeugt und aufrechterhalten werden. Wir betrachten die grundlegenden Mechanismen für die drei Haupttypen von Gasentladungen:
- Unselbständige Gasentladungen (z.B. Geiger-Müller Zähler)
- Glimmentladungen  (z.B. Niederdruck-Leuchtmittel)
- Bogenentladungen  (z.B. Hochdruck-Leuchtmittel, Plasma-Oberflächenbehandlung)

In der zweiten Vorlesung betrachten wir:
- Grundlagen der Ladungsträger-Erzeugung und Vervielfachung
- Elektrischen Durchbruch
- Glimmentladungen 
- Korona-Entladungen 
- Bogen-Entladungen

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/02-Gasentladungen.pdf)

### Simulation

[xpdp1](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/xpdp1/) ([vactube.inc](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/xpdp1/inp/vactube.inp), [dark.inc](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/xpdp1/inp/dark.inp), [breakdown.inc](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/02-Gasentladungen/xpdp1/inp/breakdown.inp))


## 3. Ionisationsgrad des Plasmas

Wie stark ist das Plasma ionisiert? Wir kennen schwach ionisierte Plasmen (z.B. Gasentladungen bei kleinen Strömen) und stark ionisierte Plasmen (z.B. Sterninneres oder Fusionsplasmen). Bei hoher Rate der Ionisations- und Rekombinationsprozesse kann man die Besetzung der verschiedenen Teilchenzustände oft wie in einem thermischen Gleichgewicht gemäß den Gesetzen der statistischen Physik durch eine gemeinsame Temperatur beschreiben, obwohl ein Plasma eigentlich fast immer kein abgeschlossenes thermodynamisches System darstellt. Bei kleinen Wechselwirkungsraten genügt es oft, die jeweils stärksten Ionisations- und Rekombinationsprozesse zu bilanzieren (Ratengleichgewicht).

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/03-Ionisationsgrad/03-Ionisationsgrad.pdf)

## 4. Coulomb-Stoesse

In vielen Plasmen dominiert die Wechselwirkung durch die langreichweitige Coulomb-Wechselwirkung über (kurzreichweitige) atomare Anregungs- oder Ionisationsstöße. Coulomb-Stöße sind elastisch, d.h. erhalten Gesamtenergie und Gesamtimpuls der Stoßpartner, da keine inneren atomaren Anregungen erfolgen. Die Abbremsung von Ladungsträgern, Ursache für den elektrischen Widerstand des Plasmas,  erfolgt durch Ablenkung der Bewegungsrichtung. Aufgrund der langen Reichweite des elektrostatischen Potenzials dominieren Stöße mit kleinen Ablenkwinkeln. Für ein reines Coulomb-Potenzial divergiert der totale Wirkungsquerschnitt (wie aus der klassischen Mechanik bekannt) und die Abbremskraft. Für ein abgeschirmtes Potenzial wie es im quasineutralen Plasma durchweg vorliegt, können wir in üblicher Näherung die Kraftwirkung außerhalb einer Abschirmlänge (Debye-Länge) vernachlässigen und so eine effektive Stossfrequenz und damit den elektrischen Widerstand des Plasmas berechnen.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/04-Coulombstoesse/04-Coulombstoesse.pdf)

### Simulation

[coulomb.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/04-Coulombstoesse/coulomb.py), 
[coulomb_anim.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/04-Coulombstoesse/coulomb_anim1.py),
[coulomb_stat.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/04-Coulombstoesse/coulomb_stat.py)

## 5. Einzelteilchen-Bewegung

Bisher haben wir die Bewegung von Testteilchen im elektrostatischen Feld betrachtet, speziell im Zentralpotenzial (Coulomb-Stöße). Wir fügen nun ein Magnetfeld hinzu. Im homogenen Magnetfeld führt die Lorentz-Kraft zur Kreisbewegung von geladenen Teilchen (Gyration, Zyklotronbewegung). In Kombination mit einem elektrischen Feld oder einer Inhomogenität des Magnetfelds ist der Radius dieser Kreisbewegung ("Gyroradius" bzw. "Larmor-Radius") moduliert, was zu einem Bahnversatz bei jedem Umlauf führt, der sog. Teilchendrift. Wir besprechen die Ursachen für Teilchendriften. Mit Magnetfeldern lassen sich außerdem geladene Teilchen räumlich einschließen, z.B. im sogenannten "magnetischen Spiegel" (oft auch "magnetische Falle" genannt). Der Austausch kinetischer Energie in der Gyrationsbewegung senkrecht zum Magnetfeld und der Bahnbewegung parallel zum Magnetfeld führt dabei am Ort des höchsten erreichbaren Magnetfelds zur Bahnumkehr. Damit lassen sich Teilchen um eine räumliches Magnetfeldminimum herum einschließen, wenn ihre Geschwindigkeit parallel zu <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{B}}"> nicht zu hoch wird. Ein Beispiel für eine solche Spiegelkonfiguration ist das Erdmagnetfeld, welches in den sogenannten "van Allen"-Gürteln geladene Teilchen aus dem Sonnenwind speichern kann.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/05-Einzelteilchen/05-Einzelteilchen.pdf)

### Simulation

[particle.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/05-Einzelteilchen/particle.py),
[gyro_anim.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/05-Einzelteilchen/gyro_anim.py), 
[drift_anim.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/05-Einzelteilchen/drift_anim.py), 
[mirror_anim.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/05-Einzelteilchen/mirror.anim.py)

## 6. Hochfrequenz-Plasmaentladungen

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/06-HFPlasmen/06-HFPlasmen.pdf)

### Simulation

[RF_particle_anim.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/06-HFPlasmen/RF_particle_anim.py)

## 7. Kinetische Beschreibung von Plasmen

Wir beginnen nun damit, die Rückwirkung der geladenen Teilchen auf die (elektrischen und magnetischen) Felder zu beschreiben. In der Einzelteilchenbeschreibung ist dies zu komplex, da Gleichungssysteme der Ordnung <img src="https://render.githubusercontent.com/render/math?math={\color{white}N^2}"> (für Teilchenzahl N) zu lösen sind, was auch heute noch unpraktisch ist.

Die Teilchengesamtheit (klassische Teilchen) wird stattdessen als kontinuierliche Phasenraumdichte <img src="https://render.githubusercontent.com/render/math?math={\color{white}f(t,x,v)}"> ausgedrückt. Für langreichweitige Felder (wie z.B. das elektromagnetische) leiten wir eine Bewegungsgleichung her, die Vlasov-Gleichung. Sie hat die Form einer Kontinuitätsgleichung für die Phasenraumdichte, was der Erhaltung der Teilchenzahl in einem Phasenraumvolumen entspricht. Kurzreichweitige Stöße können dann immer noch (mehr oder weniger ad-hoc) als Quell- und Senkenterme hinzugeführt werden, was zur Boltzmann-Gleichung der statistischen Physik führt.

Als einfaches Beispiel für die Anwendung der Vlavov-Gleichung betrachten wir hochfrequente Elektronen-Wellen in einer Raum- und einer Geschwindigkeits-Dimension, und rein elektrostatisch. Die Coulomb-Kraft fungiert als Rückstellkraft und ist daher (anti-) parallel zur Auslenkung, es handelt sich also um longitudinale Wellen. Da diese sind Vakuum nicht ausbreitungsfähig, wohl aber im Plasma. Das heisst, solche longitudialen Wellen können nicht von ausserhalb des Plasmas angeregt werden, sondern nur durch Vorgänge im Plasma, speziell Abweichungen der Verteilungsfunktion von einer homogenen Maxwell-Verteilung. Wir berechnen die Dispersionsrelation, die im Grenzfall langer Wellen (kleiner Wellenvektor <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{k}}">) die Plasmafrequenz als Lösung hat und bei höherem <img src="https://render.githubusercontent.com/render/math?math={\color{white}k}">, abhängig vom mittleren Geschwindigkeitsquadrat, ansteigt. Dieses Resultat vergleichen wir mit einer selbstgemachten Particle-In-Cell-Rechnung, die durch das beigefügte Python-Programm realisiert wird. In der Rechnung wird für eine Anzahl diskreter Teilchen in aufeinanderfolgenden Zeitschritten jeweils ein gemeinsames Feld berechnet und dann die Teilchen in diesem Feld propagiert.

[Verlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/07-KinetischeBeschreibung/07-KinetischeBeschreibung.pdf)

### Simulation
[pic_es1.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/07-KinetischeBeschreibung/pic_es1.py) 

[Blog](https://medium.com/swlh/create-your-own-plasma-pic-simulation-with-python-39145c66578b)

## 8. Landau-Dämpfung

Die Lösung der Vlasov-Gleichung erfolgte in der vorherigen Stunde durch einen Ansatz ebener Wellen, ~ <img src="https://render.githubusercontent.com/render/math?math={\color{white}e^{(i\vec{k}\cdot\vec{x} - i\omega)}}">, unendlich ausgedehnt und mit unendlicher Lebensdauer. Dieser Ansatz führt jedoch beim Aufstellen der Dispersionsrelation <img src="https://render.githubusercontent.com/render/math?math={\color{white}\omega(k)}"> zu einer Singularität für diejenigen Teilchen in der Verteilung, deren Geschwindigkeit <img src="https://render.githubusercontent.com/render/math?math={\color{white}v}"> gleich  der Phasengeschwindigkeit <img src="https://render.githubusercontent.com/render/math?math={\color{white}v_{ph} = \omega/k}"> der Welle ist (Resonanz zwischen Teilchen und Welle). Wenn solche Teilchen in der Verteilungsfunktion nicht vernachlässigt werden können, dann kann die ansonsten so beliebte ebene Welle nicht angesetzt werden, und es muss ein anderer Ansatz herhalten. Wir besprechen die klassische Lösung von Lev Landau, der den Ansatz um die Möglichkeit des Anwachsens bzw. der Dämpfung der Welle erweitert, ~ <img src="https://render.githubusercontent.com/render/math?math={\color{white}e^{(i\vec{k}\cdot\vec{x} - \gamma - i\omega)}}">, wobei gamma die sogenannte Anwachsrate ist. Die Lösung erfolgt dann mit Hilfe der Laplace-Transformation anstelle der Fourier-Transformation. Die Laplace-Transformierte ist i.a. komplex, was die Methoden der komplexen Analysis für die Lösung eröffnet, speziell die Vereinfachung der zur Rück-Transformation notwendigen Pfad-Integrale im komplexen (gamma, omega)-Raum gerade durch die Existenz von Singularitäten.

Als Beispiel untersuchen wir wiederum den Fall der elektrostatischen jeweils in Raum und Geschwindigkeit eindimensionalen Welle. Durch die Anwesenheit der Resonanz ergibt sich für monoton in v abfallende Geschwindigkeitsverteilungen eine Dämpfung, ganz ohne zusätzliche Stoßterme in der kinetischen Gleichung. Diese Dämpfung kommt durch Netto-Energieübertrag aus dem elektrischen Feld in die Teilchenverteilung um die Resonanz herum zustande, welche die Verteilungsfunktion für im elektrischen Potenzialtrog der Welle gefangene Teilchen lokal abflacht. Die Wellen-Amplitude wird  aber nicht beliebig klein - die Abflachung der Verteilungsfunktion verringert (nicht-linear) die Dämpfungsrate und eine endliche Feldamplitude bleibt übrig. Wir benutzen das Particle-In-Cell-Programm der vorigen Stunde um dieses Verhalten in einer (einfachen) numerischen Rechnung zu bestätigen.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/08-LandauDaempfung/08-LandauDaempfung.pdf)

### Simulation
[landaudamp.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/08-LandauDaempfung/landaudamp.py)

## 9. Flüssigkeitsbeschreibung von Plasma

Wir haben die kinetische Beschreibung von Plasmen mit einer Phasenraumdichte im Orts- und Geschwindigkeitsraum benutzt, um Effekte zu beschreiben, die die Geschwindigkeitsverteilung beeinflussen, speziell die Dämpfung von elektrostatischen Wellen. Für viele Phänomene jedoch muss die Geschwindigkeitsverteilung nicht vollständig berechnet werden, sondern kann durch (wenige) Momente, Integrale über <img src="https://render.githubusercontent.com/render/math?math={\color{white}(v^n \cdot f) d^3v, n=0, 1, 2, ...,}"> genügend genau wiedergegeben werden. Das führt zu wenigen, abzählbaren, Bewegungsgleichungen für Momente der Boltzmann-Gleichung, die die Form von Kontinuitätsgleichungen haben: Teilchentransport (<img src="https://render.githubusercontent.com/render/math?math={\color{white}n=0}">), Impuls/Kraftgleichung (<img src="https://render.githubusercontent.com/render/math?math={\color{white}n=1}">), Wärmetransport  (<img src="https://render.githubusercontent.com/render/math?math={\color{white}n=2}">) etc. Jede der Gleichungen nimmt Bezug auf das nächst höhere Moment der Verteilungsfunktion, was strenggenommen einem unendlichen Regress entspräche und die gewünschte Vereinfachung untergräbe. Hier werden, je nach Problem, von außen eingebrachte, mehr oder weniger heuristische Zusatzbedingungen für das abgebrochene Moment verwendet, die sogenannte Schliessung. Der so eingeführte Formalismus gilt zunächst getrennt für jede Teilchensorte. Als Anwendung werden wir in der nächsten Vorlesung elektromagnetische Mehr-Flüssigkeits-Wellen behandeln.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/09-Fluidbeschreibung/09-Fluidbeschreibung.pdf)

## 10. Wellen im kalten Plasma

Die in der vorigen Vorlesung eingeführte Flüssigkeitsbeschreibung des Plasmas vereinfacht Fälle, in denen keine subtile Störung der Geschwindigkeit-Verteilungsfunktion behandelt werden muss. Als Beispiel betrachten wir nun elektromagnetische Wellen im kalten Plasma. Die Behandlung beinhaltet ein konstantes magnetisches "Hintergrund"-Feld, elektromagnetische Wellen beliebiger Polarisationsrichtung und wahlweise Elektronen-Schwingungen sowie kombinierte Elektronen- und Ionen-Schwingungen. Einzig die Temperatur, und damit der Plasmadruck sind als vernachlässigbar angenommen, was die Behandlung mit ausschliesslich der Kraftgleichung ermöglicht. Diese Vereinfachung hat lediglich den Nachteil, dass die Dämpfung von Wellen an Resonanzstellen nicht beschrieben werden kann.

Wir behandeln die Dispersionsrelation <img src="https://render.githubusercontent.com/render/math?math={\color{white}\omega(k)}">, bzw. Phasengeschwindigkeit <img src="https://render.githubusercontent.com/render/math?math={\color{white}v_{ph}=\omega/k}"> und Brechungsindex <img src="https://render.githubusercontent.com/render/math?math={\color{white}N= c/v_{ph}}"> für die Ausbreitung parallel zum statischen Hintergrund-Magnetfeld für die zwei Zweige parallel <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{B}}">, die "R"-Welle und die "L"-Welle. Mit Hilfe dieser Ergebnisse lassen sich der "Whistler"-Mode in der Ionosphäre sowie die Faraday-Rotation der Polarisierung erklären. Weiterhin wird die Ausbreitungsrichtung senkrecht zum "Hintergrund"-Magnetfeld im kalten Plasma behandelt. Je nach Polarisationsrichtung <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{E_1}}"> parallel zu <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{B_0}}"> ("O"-Welle)  und <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{E_2}}"> senkrecht zu <img src="https://render.githubusercontent.com/render/math?math={\color{white}\vec{B_0}}"> ("X"-Welle) ergeben sich unterschiedliche Dispersionsrelationen.

[Vorleung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/10-Wellen/10-Wellen.pdf)

### Simulation
[cold_plasma_waves.py](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/10-Wellen/cold_plasma_waves.py)

## 11. Magnetohydrodynamik

Die Flüssigkeitsbeschreibung des Plasmas (ein diskreter Satz von Momentengleichungen je Teilchenspezies) kann weiter vereinfacht werden, indem diese zu einem gemeinsamen Satz Gleichungen für alle Spezies zusammengefasst werden (Einflüssigkeits-Modell, bzw. Magnetohydrodynamik, MHD). Neben der mittleren Massenströmungsgeschwindigkeit tritt nun (aufgrund der Differenzgeschwindigkeit von Ionen und Elektronen) die elektrische Stromdichte auf. An die Stelle zweier Kraftgleichungen für Elektronen und Ionen treten nun eine Kraftgleichung und das sog. verallgemeinerte Ohm'sche Gesetz, die beide jeweils die Massenströmungsgeschwindigkeit und die elektrische Stromstärke enthalten.

Mit der MHD-Beschreibung können wir nun die vorherige Vereinfachung des kalten Plasmas aufheben und die Ausbreitung von Wellen in kompressiblen Plasmen (endliche Rückstellkraft durch Druckänderung) behandeln. Die Rückstellkraft der Plasmaschwingung wird bei Auslenkung parallel B (longitudinale Welle) durch die Druckänderung (akustische Welle) und bei Auslenkung senkrecht B (transversale Welle) durch die Magnetfeldspannung (Alfvén-Welle) bewirkt.

[Vorlesung](/literature/Suttrop%20-%20Einfuehrung%20in%20Plasmaphysik/11-Magnetohydrodynamik/11-Magnetohydrodynamik.pdf)
