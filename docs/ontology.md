# Food ontology for RDD Metabolomics

## Overview

The Global FoodOmics food ontology provides a hierarchical classification system for organizing food reference data used in RDD metabolomics analyses. This ontology enables queries at varying levels of specificity, from broad food origin types (plant/animal/fungi/mineral/algae) down to specific food items and varieties.

## Hierarchical structure

The food ontology consists of **6 hierarchical groups** (sample_type_group1 through sample_type_group6), where each group provides increasingly specific classification:

```
Group 1: Origin type (plant/animal/fungi/mineral/algae)
    ↓
Group 2: Major food category
    ↓
Group 3: Food type/class
    ↓
Group 4: Food subtype
    ↓
Group 5: Specific food item
    ↓
Group 6: Specific variety or item (repeat Group 5 if not more specific)
```

## Group definitions

### Group 1: Origin type
The broadest classification based on biological or material origin.

**Categories:**
- `plant` - Plant-based foods
- `animal` - Animal-based foods
- `fungi` - Fungal-based foods (mushrooms, yeast)
- `mineral` - Mineral-based foods (salt, sparkling water)
- `algae` - Algae-based foods (seaweed, spirulina)

### Group 2: Major food category
Primary food groups organized by origin type. Detailed hierarchies are provided in the subgroup hierarchy reference.

**Plant-based categories:**
- `fruit` - All fruit types
- `vegetable/herb` - Vegetables, herbs, spices
- `sap` - Tree saps and syrups

**Animal-based categories:**
- `white meat` - Poultry and seafood
- `red meat` - Beef, pork, lamb, mutton, horse
- `game meat` - Wild game (deer, wild boar, hare, geese)
- `egg` - Eggs from various sources
- `honey` - Honey and bee products
- `dairy` - Milk, yogurt, cream, cheese
- `insect` - Edible insects
- `gastropod` - Snails and related
- `reptile` - Reptile meat

**Fungi-based categories:**
- `fungi` - Mushrooms and yeast

**Mineral-based categories:**
- `mineral` - Minerals and water

### Group 3: Food type/class
More specific food categories within major groups.

**Examples:**

**Within fruit:**
- `fleshy fruit` - Fruits with fleshy edible parts
- `achene` - Dry fruits (sunflower seeds)
- `siliqua` - Pod-like fruits (radish)
- `legume` - Legumes and beans
- `nut` - Tree nuts
- `(caryopisis)` - Grains (technical botanical term)
- `grain/grass` - Grasses and grain crops

**Within vegetable/herb:**
- `vegetable/herb` - General vegetables
- `herb/spice` - Culinary herbs and spices
- `tea/infusion` - Tea leaves and infusion materials
- `flower` - Edible flowers

**Within white meat:**
- `poultry` - Chicken, duck, turkey
- `seafood` - Marine animals (octopus, calamari/squid)
- `shell fish` - Shellfish
- `fish` - Fish species

**Within egg:**
- `egg_fish` - Fish eggs/roe
- `egg_poultry` - Bird eggs

**Within dairy:**
- `milk, yogurt, cream` - Dairy liquids and semi-solids
- `cheese` - Cheese products

**Within fungi:**
- `mushroom, yeast` - Fungal foods

**Within mineral:**
- `mineral` - Mineral products

### Group 4: Food subtype
Specific food subtypes or botanical classifications.

**Examples:**

**Within fleshy fruit:**
- `berry` - True berries (tomato, grapes, eggplant, banana, avocado, kiwi, coffee, currant, passionfruit, pepper, blueberry, guava, papaya, cranberries)
- `pepo` - Melon family (watermelon, cucumber, squash, cantaloupe, pumpkin, spaghetti squash, honeydew)
- `drupe` - Stone fruits (apricots, olives, peaches, plums, cherries, mangos, amlas, coconut, date, sapote, nectarine)
- `drupe_aggregate` - Aggregate stone fruits (blackberry, raspberry, loganberry, cherimoya, rose hip, sugar apple, custard apple, soursop)
- `pome` - Pome fruits (apple, pear, quince, loquat, medlar, ca holly, serviceberry, chokeberries, hawthorns, rowan, asian pear)
- `(hesperidium)` - Citrus 
- `citrus` - Citrus fruits (orange, citron, grapefruit, kumquat, lemon, lime)
- `pseudocarp` - False fruits (strawberries)
- `(multiflower, multiple fruits)` - Compound fruits
- `multifruit` - Multiple fruits (pineapple, fig, mulberry, breadfruit, jackfruit)

**Within grain/grass:**
- `grain/grass` - Grain types
- `wheat` - Wheat varieties
- `rice` - Rice varieties
- `barley` - Barley
- `rye` - Rye
- `hay` - Hay/fodder grasses
- `corn` - Corn/maize

**Within vegetable/herb:**
- `vegetable` - General vegetables (beet, potato, kale, arugula, spinach, lettuce, chard, broccoli, carrot)
- `herb/spice` - Herbs and spices (cilantro, mint, paprika, etc.)
- `tea/infusion` - Teas (tea, infusion)
- `flower` - Edible flowers

**Within achene:**
- `achene` - Sunflower seed

**Within siliqua:**
- `siliqua` - Radish

**Within legume:**
- `legume` - Beans (tamarind, bean)

**Within nut:**
- `nut` - Tree nuts (chestnut, hazelnut, acorn, nutmeg, walnut, almond, pecan)

**Within white meat:**
- `poultry` - Poultry types
- `seafood` - Marine animals (octopus, calamari/squid)
- `shell fish` - Shellfish (shrimp, scallop, crab, clams, mussels, winkles)
- `fish_saltwater` - Saltwater fish (shark, snapper, polluck, mackerel, herring, flounder, smelt, salmon, sardine)
- `fish_freshwater` - Freshwater fish (tilapia, trout)

**Within red meat:**
- `lamb, cow, pig, mutton, horse` - Red meat animals
- `pork, beef, chops` - Red meat cuts

**Within game meat:**
- `deer, wild boar, hare, geese` - Large game
- `small birds, rabbit` - Small game

**Within egg:**
- `egg_salmon` - Salmon roe
- `egg_chicken, egg_duck` - Poultry eggs
- `egg_poultry` - General poultry eggs
- `egg_fish` - Fish roe

**Within honey:**
- `honey` - Honey varieties

**Within dairy:**
- `milk, yogurt, cream` - Dairy liquids
- `cheese` - Cheese types

**Within insect:**
- `worm` - Edible worms

**Within gastropod:**
- `snail` - Edible snails

**Within reptile:**
- `snake, lizard` - Reptile types

**Within fungi:**
- `mushroom, yeast` - Fungal varieties

**Within mineral:**
- `mineral` - Mineral types

### Group 5: Specific food item
Individual food items or specific varieties.

**Examples:**

**Fruits:**
- Berry group: `tomato`, `grapes`, `eggplant`, `banana`, `avocado`, `kiwi`, `coffee`, `currant`, `passionfruit`, `pepper`, `blueberry`, `guava`, `papaya`, `cranberries`
- Pepo group: `watermelon`, `cucumber`, `squash`, `cantaloupe`, `pumpkin`, `spaghetti squash`, `honeydew`
- Drupe group: `apricots`, `olives`, `peaches`, `plums`, `cherries`, `mangos`, `amlas`, `coconut`, `date`, `sapote`, `nectarine`
- Drupe aggregate: `blackberry`, `raspberry`, `loganberry`, `cherimoya`, `rose hip`, `sugar apple`, `custard apple`, `soursop`
- Pome group: `apple`, `pear`, `quince`, `loquat`, `medlar`, `ca holly`, `serviceberry`, `chokeberries`, `hawthorns`, `rowan`, `asian pear`
- Citrus: `orange`, `citron`, `grapefruit`, `kumquat`, `lemon`, `lime`
- Pseudocarp: `strawberries`
- Multifruit: `pineapple`, `fig`, `mulberry`, `breadfruit`, `jackfruit`

**Seeds and grains:**
- `sunflower seed`
- `farro`, `spelt` (wheat varieties)
- `jasmine rice`, `basmati rice`, `brown rice` (rice varieties)

**Legumes:**
- `soy bean`, `green beans`, `mung bean`, `tamarind`, `bean`

**Nuts:**
- `chestnut`, `hazelnut`, `acorn`, `nutmeg`, `walnut`, `almond`, `pecan`

**Vegetables:**
- `beet`, `potato`, `kale`, `arugula`, `spinach`, `lettuce`, `chard`, `broccoli`, `carrot`, `radish`

**Herbs/spices:**
- `cilantro`, `mint`, `paprika`, `chocolate mint` (specific variety)

**Teas:**
- `green tea`, `black tea`

**Flowers:**
- `nasturtium`

**Poultry:**
- `chicken`, `duck`, `turkey`

**Seafood:**
- `octopus`, `calamari/squid`, `shrimp`, `scallop`, `crab`, `clams`, `mussels`, `winkles`

**Fish:**
- Saltwater: `shark`, `snapper`, `polluck`, `mackerel`, `herring`, `flounder`, `smelt`, `salmon`, `sardine`
- Freshwater: `tilapia`, `trout`

**Red meat:**
- `pork`, `beef`, `chops`, `lamb`, `cow`, `pig`, `mutton`, `horse`

**Game:**
- `deer`, `wild boar`, `hare`, `geese`, `small birds`, `rabbit`

**Eggs:**
- `salmon calviar`, `egg_chicken`, `egg_duck`, `egg_quail`, `egg_trout`

**Honey:**
- `honey_wildflower`, `honey_avocado`

**Dairy:**
- `milk_cow`, `milk_goat`, `yogurt_cow`, `cheese/or blue`

**Insects:**
- `worm_nightcrawler`

**Gastropods:**
- `snail`

**Reptiles:**
- `snake_rattle`, `lizard`

**Fungi:**
- `king oyster mushroom`, `crimini mushroom`, `shitake mushroom`, `maitake mushroom`, `yeast`

**Minerals:**
- `mineral`

### Group 6: Specific variety or item
The most specific level.

**Rule:** If sample is not more specific than Group 5, repeat the value.

**Example:**
```
Group 5: strawberries
Group 6: strawberries  (repeated because no more specific variety)

Group 5: basil
Group 6: purple basil
```


## Important notes

## Metadata file format

The food ontology is encoded in a .txt (separated by tabs) file with the following structure:

```tsv
filename	sample_type_group1	sample_type_group2	sample_type_group3	sample_type_group4	sample_type_group5	sample_type_group6
apple_red_001.mzML	plant	fruit	fleshy fruit	pome	apple	apple
orange_002.mzML	plant	fruit	fleshy fruit	citrus	orange	orange
salmon_003.mzML	animal	white meat	fish	fish_saltwater	salmon	salmon
strawberry_004.mzML	plant	fruit	fleshy fruit	pseudocarp	strawberries	strawberries
```

**Column requirements:**
- `filename`: Must match the reference spectrum filename in GNPS
- `sample_type_group1` through `sample_type_group6`: Hierarchical classifications following the ontology structure
- **Group 6 rule**: If no more specific than Group 5, repeat the Group 5 value
- Use consistent naming as defined in the ontology (case-sensitive)

## Usage in RDD analysis

### Understanding groups vs. levels

In the RDD package, the `level` parameter maps to groups as follows(Using GFOP ontology):
- `level=0` → Individual files (no grouping by ontology)
- `level=1` → Group 1 (sample_type_group1: Origin type)
- `level=2` → Group 2 (sample_type_group2: Major food category)
- `level=3` → Group 3 (sample_type_group3: Food type/class)
- `level=4` → Group 4 (sample_type_group4: Food subtype)
- `level=5` → Group 5 (sample_type_group5: Specific food item)
- `level=6` → Group 6 (sample_type_group6: Specific variety)

### Selecting analysis levels

Choose the appropriate level based on your research question:

**Individual file level (level 0)**
```python
# Analyze each reference file individually
rdd.get_counts(level=0)
```

**Broad dietary patterns (levels 1-2)**
```python
# Compare plant vs. animal food consumption (Group 1)
rdd.get_counts(level=1)

# Analyze major food groups: fruits, vegetables, meats, dairy (Group 2)
rdd.get_counts(level=2)
```

**Intermediate classification (levels 3-4)**
```python
# Examine specific food types: fleshy fruits, grains, poultry (Group 3)
rdd.get_counts(level=3)

# Detailed food subtype analysis: berries, citrus, fish types (Group 4)
rdd.get_counts(level=4)
```

**Fine-grained analysis (levels 5-6)**
```python
# Specific food items: strawberries, salmon, wheat (Group 5)
rdd.get_counts(level=5)

# Most specific varieties and items (Group 6)
rdd.get_counts(level=6)
```

### Multi-level analysis

Compare patterns across multiple levels to understand dietary composition at different granularities:

```python
# Broad overview
broad = rdd.get_counts(level=2)  # Major food categories

# Detailed breakdown
detailed = rdd.get_counts(level=4)  # Food subtypes

# Visualize hierarchical flow
viz.plot_sankey_diagram(rdd, start_level=1, end_level=4)
```

## Extending the ontology

RDD metabolomics works best when your reference data is organized with a proper hierarchical ontology. The structure and number of groups can be adapted to fit your specific research domain and the natural classification system of your reference materials. 

When creating custom reference metadata:
1. Define clear hierarchical groups that make sense for your domain
2. Maintain consistent naming conventions across all reference samples
3. Ensure each lower group nests properly within higher groups
4. Apply the repeat rule for the final group when appropriate
5. Document your ontology structure for reproducibility

## Additional resources

- **Full dataset**: MassIVE repository MSV000084900
- **GNPS tutorial**: https://ccms-ucsd.github.io/GNPSDocumentation/tutorials/rdd/
- **Package documentation**: https://gnps-rdd.readthedocs.io
