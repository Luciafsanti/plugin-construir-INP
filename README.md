# plugin-construir-INP — Plugin de QGIS

Genera un modelo **EPA SWMM 5.2** de escurrimiento superficial en calles (sistema mayor, conductos con sección `STREET`) a partir de:

- **DEM** (raster)
- **Polígono de área de interés / extent**
- **Ejes de calle** (vector de líneas)
- **Manzanas** (vector de polígonos)
- **Uso de suelo clasificado** (raster)
- **Hietograma / serie de lluvia** (CSV)

Produce:
- Un archivo **`.inp`** listo para abrir en SWMM 5.2.
- Tres capas de salida en QGIS (Nodos, Conductos, Subcuencas) — la representación GIS equivalente al `.inp`.

Implementa como plugin de Processing las recomendaciones del workflow: nodos en intersecciones/puntos bajos, conductos de calle con sección `STREET`, subcuencas delimitadas por manzana y asignadas al nodo más cercano, corrección de manzanas grandes, e infiltración tipo Horton por defecto.

---

## 1. Instalación

1. Abrí QGIS → **Complementos → Administrar e instalar complementos → Instalar desde ZIP**.
2. Seleccioná el archivo `swmm_street_builder.zip`.
3. Activá el plugin cuando lo pida el instalador.
4. Abrí la **Caja de herramientas de Procesos** (Processing Toolbox) → buscá **"SWMM Street Builder"** → algoritmo **"Construir modelo SWMM de escurrimiento en calles"**.

> Alternativa manual: copiar la carpeta `swmm_street_builder/` completa dentro de tu carpeta de plugins de QGIS (`~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/` en Linux, o el equivalente en Windows/Mac) y reiniciar QGIS.

---

## 2. Antes de correrlo

### 2.1 Preparación de capas
- Todas las capas de entrada deben estar en el **mismo CRS proyectado** (metros).
- La capa de ejes de calle **no necesita venir pre-noded**: el plugin la parte automáticamente en cada cruce y agrupa (clusteriza) los extremos cercanos como nodos.
- El raster de uso de suelo debe ser **clasificado** (valores enteros de categoría, no continuos).

### 2.2 Editar las tablas de referencia (`lookup_tables/`)
- **`street_categories.csv`**: define los parámetros de la sección `STREET` (Tcrown, Hcurb, Sx, nRoad, etc.) por categoría de calle. Vienen 3 categorías de ejemplo (`local`, `colectora`, `arterial`) con valores típicos — **ajustalos a tu normativa vial local**. Si tu capa de ejes tiene un campo de categoría, indicalo en el parámetro "Campo de categoría de calle"; si no lo indicás, se aplica la primera categoría de la tabla a todos los tramos.
- **`landuse_imperviousness.csv`**: relaciona cada código de clase de tu raster de uso de suelo con un % de impermeabilidad (`pct_imperv`) y con un **Número de Curva SCS (`cn`)**. La columna `cn` la usa solo el módulo de drenaje natural (sección 3.4) para calcular la lluvia efectiva. **Los valores de ejemplo (CN de grupo hidrológico B, condición AMC II) son solo para probar el plugin — reemplazalos por los códigos y CN reales de tu clasificación.** Si falta la columna `cn` se usa CN 75 con un warning.
- **`hietograma_ejemplo.csv`**: formato de ejemplo para el parámetro de lluvia (columnas `time,value`, o `date,time,value` si tu serie cubre más de un día). `value` se interpreta como intensidad (mm/hr) por defecto.

---

## 3. Parámetros del algoritmo

| Parámetro | Qué hace |
|---|---|
| DEM / Extent / Ejes de calle / Manzanas / Uso de suelo | Las 5 capas de entrada obligatorias |
| Campo de categoría de calle | Campo (opcional) en Ejes de calle con el nombre de categoría (debe matchear la columna `category` del CSV) |
| Tabla de categorías de calle (CSV) | Ver 2.2 |
| Tabla uso de suelo → %impermeable (CSV) | Ver 2.2 |
| Hietograma por defecto (CSV) | Serie de lluvia única, se usa **solo si no se cargan pluviómetros** (ver 3.1) |
| Pluviómetros / Rain Gauges (shp puntos, opcional) | Ver 3.1 |
| Campo de código de pluviómetro | Campo en la capa de pluviómetros que matchea la columna `codigo` del CSV |
| Tabla de pluviómetros (CSV, opcional) | Ver 3.1 |
| Ríos / cauces receptores (shp lineal, opcional) | Ver 3.2 — son **salidas** del modelo, no conductos |
| Embalses / lagunas receptores (shp polígonos, opcional) | Ver 3.3 — son **salidas** del modelo, no Storage Units |
| Tolerancia de unión de nodos (m) | Distancia bajo la cual dos extremos de calle/río se consideran el mismo nodo. Con DEM de baja resolución o capas mal snapeadas, subilo (1–2 m) |
| Área máxima de subcuenca antes de subdividir (ha) | Manzanas más grandes que este umbral se subdividen con una grilla regular |
| Coeficiente Width = Area/(coef·Perímetro) | Aproximación clásica de ancho de subcuenca cuando no hay longitud de escurrimiento medida |
| Distancia máx. calle→río/storage para crear outfall (m) (def. 20) | Un nodo de **grado 1** (extremo de calle) a ≤ esta distancia de un río o storage se marca como `OUTFALL`. Ya no se usa el borde del AOI |
| Título del modelo | Texto libre para `[TITLE]` |

### 3.1 Pluviómetros múltiples (opcional)

Si tu ciudad tiene más de un pluviómetro/hietograma, cargá:

- **Capa de puntos** (uno por pluviómetro), con un campo de código (el que indiques en "Campo de código de pluviómetro").
- **CSV de asociación**, con columnas exactas:

| codigo | rain_format | file |
|---|---|---|
| PLUV01 | intensity | hietograma_pluv01.csv |
| PLUV02 | cumulative | hietograma_pluv02.csv |

- `codigo`: debe matchear el valor del campo de código en la capa de puntos.
- `rain_format`: `intensity`, `volume` o `cumulative` (minúsculas, tal como indicaste) — se traduce internamente a `INTENSITY`/`VOLUME`/`CUMULATIVE` de SWMM.
- `file`: ruta al CSV del hietograma de ese pluviómetro (mismo formato `time,value` o `date,time,value` que el hietograma por defecto). Puede ser ruta absoluta o relativa a la ubicación del CSV de asociación.

**Ponderación de Thiessen.** La lluvia de cada subcuenca (urbana) y de cada cuenca de drenaje natural se pondera por **polígonos de Thiessen**: se calcula la fracción de área de cada elemento que cae en la zona de cada pluviómetro. Como SWMM admite un solo pluviómetro por subcuenca, la ponderación se resuelve antes de escribir el `.inp`:

- Los pesos se redondean a pasos de **5 %** (parámetro *resolución de los pesos de Thiessen*) y las estaciones que quedan en 0 se descartan.
- Si queda una sola estación, la subcuenca usa **el pluviómetro real**.
- Si quedan varias, se genera un hietograma ponderado `Σ wᵢ·Pᵢ(t)` (volumen por paso, mm) en una grilla común con el paso más fino. Las series con otro paso o inicio se reparten uniformemente dentro de cada paso, conservando el volumen. Fuera de su registro, una estación aporta 0. En el `.inp` aparece como pluviómetro sintético `RG_TH<n>` (formato `VOLUME`) con su serie `TS_TH<n>`, uno por cada combinación distinta de pesos.
- La capa de subcuencas agrega el campo `pesos` (ej. `RG_PLUV01:0.60|RG_PLUV02:0.40`). Capa opcional de salida: polígonos de Thiessen recortados al AOI.
- Las cuencas naturales usan el hietograma ponderado para el HU y calculan la lluvia efectiva (SCS-CN) **después** de promediar la lluvia.
- El intervalo de cada pluviómetro en `[RAINGAGES]` se toma de las marcas de tiempo de su serie (antes estaba fijo en 0:05).
- Los pluviómetros sin serie en la tabla CSV se omiten y su área pasa a las estaciones vecinas.

Si no cargás esta capa, se usa un único pluviómetro (`RG1`) con el "Hietograma por defecto" para todas las subcuencas, igual que antes.

### 3.2 Ríos / cauces receptores (opcional)

Cargá una capa de líneas con los ríos o canales que **reciben** el agua del modelo. **Ya no se modelan como conductos** (ni transectos): se usan solo como referencia de salida, recortados al AOI. Los campos del `.shp` (`manning_n`, `bottom_w`, `profile_csv`, etc.) se ignoran.

Dónde se crean los `OUTFALL`: solo un extremo de calle de **grado 1** que esté a ≤ "Distancia máx. calle→río/storage para crear outfall" de un río o storage. Un `OUTFALL` de SWMM admite un único conducto, por eso una calle que **cruza** un río sin terminar ahí no genera outfall. Si no hay ningún candidato, se usa como último recurso el nodo de grado 1 de menor cota (con advertencia en el log).

### 3.3 Embalses / lagunas receptores (opcional)

Cargá una capa de polígonos. Igual que los ríos, son **salidas**: no se crean `Storage Units` (`[STORAGE]` queda vacío) ni se fusionan nodos. Los polígonos se **restan a las calles** antes de nodar, de modo que las calles terminan en la orilla; los extremos de calle de grado 1 a ≤ el buffer de outfalls pasan a ser `OUTFALL`. Los campos (`ymax_m`, `y0_m`, `fevap`) ya no se usan.

### 3.4 Red de drenaje natural (opcional)

Cargá una capa de líneas (arroyos, cauces, quebradas). Módulo **opcional e independiente** (`natural_drainage.py`): si no cargás la capa, el algoritmo funciona igual que antes.

**Los cauces naturales NO son conductos SWMM ni sus cuencas son subcuencas SWMM.** La red y el DEM solo sirven para (a) decidir qué drenajes aportan a la ciudad y (b) calcular el hidrograma que entra a la red de calles como `INFLOW`.

**Qué hace:**
1. Recorta la capa al AOI y **descuenta lo que cae dentro de manzanas** (área urbana ya controlada por el drenaje de calles). Parte los ejes en cruces y confluencias y **orienta cada tramo aguas arriba → aguas abajo** según el DEM.
2. **Outlets**: el extremo aguas abajo que ningún otro tramo continúa. Solo los outlets se vinculan a la calle; los tramos intermedios **no** se cortan aunque pasen cerca de una calle.
3. **Purga y salidas**: si el outlet está a ≤ "Distancia máxima (buffer) para vincular un drenaje natural con un nodo de calle" de un **río o storage** y más cerca que el extremo de calle más cercano, el drenaje desemboca en una **salida del modelo**: se descarta (sin INFLOW) y sus tramos aguas arriba se ven en la capa de red como `tipo = a_salida`. Si no, se elimina todo drenaje cuyo outlet no tiene un nodo de calle dentro del buffer (`tipo = purgado`). El `.shp` original no se modifica.
4. Cada outlet conectado se asigna al **nodo de calle más cercano** (no OUTFALL) dentro del buffer.
5. **Cuenca de aporte** de cada outlet desde el DEM, con numpy (sin GRASS/SAGA): relleno de depresiones (priority-flood), dirección de flujo D8, el DEM se "quema" 1 m bajo los ejes conservados para que el flujo siga el cauce, y el outlet se ajusta a la celda de mayor acumulación entre las celdas de toda la red que drena a él (respaldo: radio de 2 celdas). El dominio es el AOI menos las manzanas (una divisoria que llega al borde del AOI queda recortada por el AOI y se avisa).
6. **Parámetros** de la cuenca: área, `L` = camino de flujo más largo hasta la divisoria, `H` = cota máxima de la cuenca − cota del outlet, y `CN` ponderado por área desde la cobertura (columna `cn` de `landuse_imperviousness.csv`).
7. **Tc — California Culverts Practice**: `Tc = 60·(11.9·L³/H)^0.385` con `L` en **millas** y `H` en **pies** (la conversión desde metros se hace en `tc_california_culverts()`); resultado en minutos.
8. **HU triangular SCS** con el **paso de tiempo del hietograma** ponderado por Thiessen de los pluviómetros que cubren la cuenca (o el hietograma general): `t_lag = 0.6·Tc`, `Tp = Δt/2 + t_lag`, `Tb = 2.67·Tp`, `Qp = 0.208·A[km²]/Tp[h]` en m³/s por mm (ordenadas normalizadas para conservar el volumen).
9. **Lluvia efectiva** por SCS-CN sobre la lluvia acumulada (`Ia = λ·S`, λ = parámetro *Ia/S*, 0.2 por defecto). Si la lluvia total no supera Ia la cuenca **no genera escorrentía**: no se crea INFLOW, pero la cuenca igual se dibuja en la capa de verificación (`inflow = 0`) y **convolución** `Q[n] = Σ Pe[m]·U[n−m]` (N = M + K − 1 ordenadas).
10. En el `.inp`: una serie `UH_<nodo>` en `[TIMESERIES]` con `Q / Pe_total` (m³/s por mm efectivo total) y en `[INFLOWS]`: `nodo FLOW UH_<nodo> FLOW 1.0 <Pe_total>` (Scale Factor = lluvia efectiva total en mm). Si varios outlets caen en el mismo nodo, SWMM admite un solo INFLOW `FLOW` por nodo, así que se suman los hidrogramas en una única serie con Scale Factor 1.0.

Tipos de hietograma soportados (según `rain_format` del pluviómetro): `volume` (mm por paso), `intensity` (mm/h) y `cumulative` (mm). El paso de tiempo sale de las marcas de tiempo del hietograma (no es fijo).

**Capas de salida opcionales (solo verificación, no entran al `.inp`):**
- *Cuencas de aporte*: polígonos con `node, gage, area_ha, L_m, H_m, Tc_min, CN, P_mm, Pe_mm, dt_min, Qp_unit, Qmax_cms`.
- *Red de drenaje natural*: `tipo` = `drenaje` (conservado) / `a_salida` (drena a río/storage, sin aporte) / `purgado` / `conector` (outlet → nodo) / `sin_conexion` (línea del outlet purgado hasta la calle más cercana, con su distancia en `length_m`: sirve para decidir si conviene subir el buffer).

**Importante:**
- El DEM y el raster de cobertura deben estar en el mismo CRS (métrico) que las calles.
- El hidrograma puede durar más que el `END_TIME` por defecto del `.inp` (06:00:00): revisá *Simulation Options* en SWMM (el algoritmo avisa).
- Las cuencas son de resolución del DEM (bordes en escalera): sirven para verificar divisorias, no como límites catastrales.

---

## 4. Qué hace el algoritmo, paso a paso

1. Recorta calles, manzanas, ríos, storage y pluviómetros (los que se hayan cargado) al área de interés. Los polígonos de storage se restan a las calles.
2. Limpia topología de la red de calles (fix geometries, explode multipart, auto-noding en cruces, snap). Los ríos no entran en la red.
3. Genera **nodos** agrupando extremos de calle cercanos (tolerancia configurable).
4. Muestrea la **cota de terreno** del DEM en cada nodo, promediando un patrón en cruz de 5 puntos (reduce el ruido de un DEM de resolución gruesa).
5. Detecta **outfalls**: nodos de grado 1 a ≤ el buffer de un río o storage.
6. (Sin Storage Units: los storages son salidas, `[STORAGE]` queda vacío.)
7. Si no queda ningún `OUTFALL`, intenta usar el nodo de grado 1 de menor cota de toda la red; si tampoco existe (red completamente cerrada, sin extremos), **el algoritmo se detiene con error**.
8. Asigna cada conducto de calle a una **categoría** (`[STREETS]` + `[XSECTIONS]` tipo `STREET`).
9. Genera **subcuencas** a partir de las manzanas (subdividiendo las que superan el área máxima), calcula `%Imperv` y `%Slope` por zona, asigna cada una al **nodo más cercano** como `Outlet` y al **pluviómetro más cercano** como `Raingage`.
10. Arma `[RAINGAGES]`/`[TIMESERIES]` (uno o varios, según se hayan cargado pluviómetros), y escribe el `.inp` completo más las 3 capas de salida (Nodos —incluye tipo `STORAGE`—, Conductos —incluye tipo calle/río—, Subcuencas).

---

## 5. Después de correrlo — checklist obligatorio

Este plugin automatiza el 80% del trabajo geométrico, pero **el modelo generado no debe simularse sin revisión**:

- [ ] Abrir el `.inp` en SWMM y revisar visualmente la geometría de la red.
- [ ] Verificar los **outfalls** (extremos de calle junto a ríos/storages) y que cada calle que debe descargar termine en la orilla; revisar el log (cantidad de outfalls y de outlets de drenaje que caen a salida).
- [ ] Revisar tramos con pendiente longitudinal muy baja o negativa (ruido del DEM).
- [ ] Ajustar `Tcrown`, `Sx`, `W` por categoría de calle según tu normativa vial real — el DEM no puede estimar esto de forma confiable.
- [ ] Confirmar que los códigos de `landuse_imperviousness.csv` coinciden con tu raster real.
- [ ] Correr una simulación corta y revisar el reporte de continuidad (`[REPORT]`) antes de una corrida completa.

Ver también la guía de errores comunes y checklist de calibración que armamos en la conversación (continuidad, inestabilidad numérica, paso de tiempo).

---

## 6. Limitaciones conocidas

- La subdivisión de manzanas grandes usa una grilla rectangular simple, no un criterio hidrológico de escurrimiento interno.
- El ancho de subcuenca (`Width`) es una aproximación geométrica (Área/Perímetro), no una longitud de escurrimiento medida.
- El modelo no incluye sumideros/bocas de tormenta (`[INLETS]`) ni una red de tuberías subterránea — es un modelo de sistema mayor puro. Si más adelante querés un modelo dual (calle + tubería), la estructura de nodos generada es un buen punto de partida para agregar la red de tuberías por separado.
- No se testeó dentro de una instalación real de QGIS (entorno de generación sin interfaz gráfica) — revisar el primer run con atención y reportar cualquier error de la API de Processing.
- Las cuencas de **drenaje natural** se delimitan con D8 sobre el DEM (numpy), quemando la red conservada; en terreno muy plano o con DEM ruidoso las divisorias pueden ser imprecisas — revisar la capa opcional de cuencas. Solo se usa lluvia efectiva SCS-CN (no hay pérdidas por humedad antecedente).
- Un drenaje natural cuyo outlet no encuentra ningún nodo de calle dentro del buffer se **purga** (no aporta al modelo). Ampliá la distancia si querés conservarlo.
- Ríos y storages son solo salidas: no se simula su nivel, transecto ni almacenamiento (el outfall es `FREE`). Un drenaje natural que desemboca en ellos no aporta caudal al modelo.
- Una calle que cruza un río sin terminar en él no genera outfall (SWMM admite un solo conducto por `OUTFALL`).
