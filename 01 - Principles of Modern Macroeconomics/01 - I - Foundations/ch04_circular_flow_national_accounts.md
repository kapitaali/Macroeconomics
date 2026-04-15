# Chapter 4 — The Circular Flow and National Income Accounting

> *"The national income of a country can be looked at in three ways: as the total of incomes earned, as the total of expenditures, or as the total of outputs."*
> — Richard Stone, Nobel Lecture, 1984

---

Chapter 3 introduced the major macroeconomic aggregates in isolation. This chapter asks a more structural question: how do these aggregates fit together? How does money flow through an economy — from firms to households as wages and profits, from households back to firms as consumption, from households and firms to the government as taxes, from the government and abroad back as spending? Understanding these flows is not merely bookkeeping. The accounting identities that govern the circular flow constrain what is macroeconomically possible and reveal connections between the fiscal deficit, private saving, and the current account that are not obvious from any single sector's perspective.

Richard Stone, who received the 1984 Nobel Prize in Economics for designing the modern System of National Accounts, spent much of his career showing that macroeconomics without consistent accounting is like physics without units: internally consistent reasoning built on an incoherent foundation. This chapter develops the framework Stone helped create and shows why it matters.

---

## 4.1 The Circular Flow of Income and Expenditure

The simplest macroeconomic model has two sectors: households and firms. Firms hire factor inputs from households — labor and capital — and pay for them as wages $W$ and profits $\Pi$. Households use this income to purchase the goods and services firms produce. Money flows in a circle: payments from firms to households as factor income, and payments from households to firms as consumer expenditure $C$.

This circular structure immediately implies a fundamental accounting identity: **aggregate income equals aggregate expenditure equals aggregate output**. Every dollar of expenditure is simultaneously someone's income — the payment for a good is the revenue that allows the firm to pay its workers and shareholders — which is simultaneously the value added in production. This triple equality is not a theory. It is a consequence of double-entry accounting applied at the economy-wide level, and it holds in every period, by definition.

Formally, in the two-sector model:

$$Y \equiv C + S_H,$$

where $Y$ is household income (wages plus profits) and $S_H$ is household saving — income not spent on consumption. From the firm's perspective, all output is either sold to households or unsold, and unsold output is treated as inventory investment $\Delta\text{Inv}$:

$$Y \equiv C + \Delta\text{Inv}.$$

Setting the two expressions equal gives $S_H \equiv \Delta\text{Inv}$: in the two-sector economy, private saving equals investment identically. This is the seed of the saving-investment identity.

### The Four-Sector Economy: Injections and Leakages

In a more realistic economy with a government sector and an external sector, the circular flow has **injections** — flows that add to demand for domestic output — and **leakages** — flows that remove purchasing power from the domestic economy.

**Injections:** investment $I$ (firms' spending on new capital goods), government purchases of goods and services $G$, and exports $X$ (foreign residents' demand for domestic output).

**Leakages:** household saving $S$ (income not spent on consumption), taxes net of transfers $T$, and imports $M$ (domestic spending on foreign rather than domestic goods).

Equilibrium in the goods market requires that injections equal leakages, which gives the standard expenditure identity:

$$Y \equiv C + I + G + (X - M) \equiv C + I + G + NX.$$

The income-side decomposition:

$$Y \equiv C + S + T + M.$$

Setting the two expressions for $Y$ equal and canceling $C$:

$$(S - I) + (T - G) = X - M \equiv NX.$$

**Definition (Saving-Investment Identity).** The **saving-investment identity** states that the sum of the private sector financial surplus $(S - I)$ and the public sector financial surplus $(T - G)$ equals net exports $NX$:

$$(S - I) + (T - G) = NX.$$

This identity holds as an accounting necessity regardless of the state of the economy and regardless of any behavioral assumption. It cannot be violated, but it can be interpreted in different directions of causation. A country running a current-account deficit ($NX < 0$) must, by arithmetic, be running a private saving deficit ($S < I$), a fiscal deficit ($T < G$), or some combination of the two. The identity does not tell us which causes which — that is a behavioral question — but it tells us that these imbalances cannot exist independently.

---

## 4.2 Sectoral Financial Balances

The saving-investment identity can be extended to a three-sector decomposition that is particularly useful for macroeconomic diagnosis and for tracking the economy's financial plumbing. Let:

- $S_H - I_H$ = household sector financial surplus (household saving minus residential investment and consumer durables)
- $S_F - I_F$ = corporate sector financial surplus (retained earnings minus business investment)
- $S_G = T - G$ = government sector financial surplus (tax revenue minus public expenditure, i.e., the fiscal surplus; negative when the government runs a deficit)
- $CA$ = current account balance (net exports plus net factor income from abroad plus net transfer payments)

The **sectoral financial balances identity**:

$$\underbrace{(S_H - I_H)}_{\text{household surplus}} + \underbrace{(S_F - I_F)}_{\text{corporate surplus}} + \underbrace{S_G}_{\text{government surplus}} = CA.$$

This identity, developed into a systematic macroeconomic framework by Wynne Godley (1999), is a powerful organizing device. Its key implication: if any one sector moves toward deficit, one or more of the remaining sectors must move toward surplus. Sectoral balances cannot all deteriorate simultaneously — the arithmetic prevents it. Conversely, if the government moves sharply toward surplus (fiscal consolidation), the private sector or the external sector must absorb the counterpart move.

### The U.S. Before the 2008 Financial Crisis: A Sectoral Balance Diagnosis

Applied to the United States in the period 2003–2007, the sectoral balance framework provides a retrospective of unusual clarity. During this period:

- The government was running a deficit: $S_G < 0$ (roughly $-3\%$ of GDP).
- Households were saving less than they were investing (including housing): $S_H - I_H < 0$ (the household sector moved to roughly $-4\%$ of GDP as mortgage borrowing surged and the personal saving rate fell toward zero).
- Because both the private and public sectors were in deficit, the current account was necessarily in deficit: $CA \approx -6\%$ of GDP.

This diagnosis does not require a behavioral model. It requires only the identity. The question of whether the current-account deficit was caused by the fiscal deficit, or by household borrowing, or by foreign capital seeking U.S. assets, is a behavioral question — but the identity tells us that all three were connected, and that any proposed explanation must be internally consistent with the accounting.

---

## 4.3 The National Income and Product Accounts

The national accounts are the systematic statistical implementation of these identities. The United States — and most countries following the SNA 2008 framework — constructs a hierarchy of aggregate income measures from GDP downward, each stripping out or adding back a specific item to move from gross production to the income actually available to households.

**Definition (Gross National Product).** **GNP** is the total income earned by a country's *residents*, regardless of where that income is generated. It differs from GDP — which measures production *within* a country's borders, regardless of the nationality of the producers — by net factor income from abroad $NFIA$:

$$\text{GNP} = \text{GDP} + \underbrace{(\text{factor income received by residents from abroad} - \text{factor income paid to foreigners domestically})}_{\text{Net Factor Income from Abroad (NFIA)}}.$$

The distinction matters for countries with large numbers of migrant workers — where earnings sent home mean GNP $>$ GDP — and for countries with heavy foreign ownership of domestic capital, where profit repatriation means GNP $<$ GDP. Ireland is a striking modern example: U.S.-headquartered multinationals book large profits in Ireland (raising Ireland's GDP), but those profits largely repatriate to the United States (so Ireland's GNP is substantially below its GDP).

**Definition (Net National Product).** **NNP** = GNP $-$ $D$, where $D$ is the consumption of fixed capital (depreciation). NNP is a better measure of *sustainable* income than GNP, because some current production is devoted to replacing worn-out capital rather than adding to productive capacity. GDP is used more commonly because depreciation is harder to measure precisely than gross investment, and because the conceptual distinction matters less for short-run fluctuations than for long-run welfare comparisons.

The complete income decomposition from GDP to personal disposable income:

$$\text{NNP} = \text{GNP} - D$$

$$\text{NI} = \text{NNP} - T^{ind} + \text{Subsidies}$$

$$\text{PI} = \text{NI} - \text{Retained earnings} - T^{corp} - T^{payroll} + \text{Transfers}$$

$$\text{PDI} = \text{PI} - T^{personal}$$

where $T^{ind}$ are indirect taxes (VAT, excise duties), $T^{corp}$ are corporate taxes, and $T^{payroll}$ are employer-side social contributions. Each step strips out income that was generated in production but does not reach households' budget constraints.

**Definition (Personal Disposable Income).** **PDI** is the income available to households after all taxes have been paid and all transfer payments received. It is the relevant budget constraint for households' consumption-saving decisions. The aggregate household saving rate $s = S_H/\text{PDI}$ is a central variable in growth theory [Ch. 5] and in debates about the effectiveness of fiscal policy [Ch. 22].

### Reconciling the Three Approaches to GDP

Chapter 3 introduced the three measurement approaches to GDP. Within the national accounts framework, their equivalence is not merely definitional but the result of explicit reconciliation procedures carried out by statistical agencies. The **expenditure approach** ($Y = C + I + G + NX$) must equal the **income approach** ($Y = W + \Pi + \text{mixed income} + T^{ind} + D$) must equal the **production approach** ($Y = \sum_j \text{VA}_j$). In practice, independent estimates of these three aggregates never agree exactly, and the residual is reported as a **statistical discrepancy** — a visible reminder that national accounts are measured with error, not derived analytically.

The income approach decomposes GDP as:

$$Y = \underbrace{W}_{\text{compensation of employees}} + \underbrace{\Pi}_{\text{gross operating surplus}} + \underbrace{MI}_{\text{mixed income (self-employment)}} + \underbrace{T^{ind} - S^{ind}}_{\text{taxes less subsidies on production}}.$$

The labor share of income, $\omega = W/Y$, is a central variable in debates about inequality [Ch. 38] and in calibrating production functions in growth models [Ch. 5]. For the United States, $\omega$ has hovered around 0.60–0.65 for most of the postwar period, with a notable downward drift since the early 2000s that has attracted substantial research attention.

---

## 4.4 The Balance of Payments

The **balance of payments (BOP)** is the system of accounts that records all economic transactions between residents of one country and residents of all other countries during a given period. Every transaction appears twice — once as a credit (something received) and once as a debit (something given), at equal values — so the overall balance of payments always sums to zero.

**Definition (Current Account).** The **current account (CA)** records: (i) the **trade balance** in goods and services (merchandise and service exports minus imports); (ii) **net primary income** — factor income earned cross-border (wages for workers, dividends and interest for capital owners); and (iii) **net secondary income** — current transfers including remittances and foreign aid. The current account balance is approximately equal to net exports $NX$ for countries with small cross-border factor income flows.

**Definition (Financial Account).** The **financial account (FA)** records net cross-border acquisitions of financial assets minus net acquisitions of financial liabilities. A positive FA (a net capital inflow, or "capital surplus") means foreigners are accumulating more claims on domestic assets than domestic residents are accumulating on foreign assets. By convention in the current SNA 2008, the financial account measures *net lending to the rest of the world*, so:

$$CA + FA = 0 \quad \text{(ignoring reserve changes and the capital account)}.$$

A current-account deficit ($CA < 0$) must be financed by a net financial inflow ($FA > 0$): foreigners are effectively lending to the domestic economy.

**Definition (Official Reserve Account).** The **official reserve account ($\Delta R$)** records changes in the central bank's holdings of foreign exchange reserves, gold, and IMF special drawing rights. An increase in reserves ($\Delta R > 0$) means the central bank is acquiring foreign assets, often to prevent currency appreciation under a managed or pegged exchange rate.

Including reserve changes and the small **capital account** (recording one-time asset transfers such as debt forgiveness), the full BOP identity is:

$$CA + FA + KA + \Delta R = 0,$$

where $KA$ is the capital account. This identity holds by construction: because every transaction is recorded twice, any failure to balance represents an **errors and omissions** term ($E\&O$) that statistical agencies report separately. Large $E\&O$ terms typically indicate unrecorded capital flows — often a sign of capital flight or offshore activity.

### The BOP and Exchange Rate Pressure

The BOP framework connects directly to exchange rate dynamics [Ch. 21, Ch. 26]. When a country runs a current-account deficit, it needs foreigners to willingly accumulate domestic financial assets. If foreigners become reluctant to do so — if the financial account surplus shrinks — the central bank must either draw down reserves or allow the currency to depreciate until the current account improves. Countries that run out of reserves without attracting sufficient capital inflows face a **sudden stop** — an abrupt reversal of capital inflows that forces a sharp, disruptive current-account adjustment. The arithmetic of the BOP makes such crises predictable in hindsight, though their timing remains notoriously difficult to forecast.

---

## 4.5 Input–Output Analysis

The production-side national accounts aggregate value added across all industries into a single GDP figure. A more disaggregated and structurally informative view is provided by the **Leontief input-output framework** (Leontief 1941), which maps the full web of inter-industry transactions.

**Definition (Technical Coefficient).** The **technical coefficient** $a_{ij}$ is the dollar value of input from industry $i$ required to produce one dollar of gross output in industry $j$. It represents the intermediate input intensity of industry $j$'s production with respect to industry $i$'s product. The $n \times n$ **technical coefficient matrix** $A = [a_{ij}]$ encodes the full structure of inter-industry dependencies.

Define the gross output vector $\mathbf{x} = (x_1, \ldots, x_n)'$ and the final demand vector $\mathbf{d} = (d_1, \ldots, d_n)'$. Final demand $\mathbf{d}$ includes private consumption, investment, government spending, and exports — the uses of output that leave the production network. The accounting identity for industry $i$: gross output equals deliveries to intermediate users plus deliveries to final demand:

$$x_i = \sum_{j=1}^n a_{ij} x_j + d_i, \quad i = 1, \ldots, n.$$

In matrix form: $\mathbf{x} = A\mathbf{x} + \mathbf{d}$, so $(I - A)\mathbf{x} = \mathbf{d}$.

**Proposition (Leontief Inverse).** If $(I - A)$ is invertible and $(I - A)^{-1} \geq 0$ element-by-element — which is guaranteed when the economy is productive, meaning each industry uses less than one dollar of inputs per dollar of output — the unique solution is:

$$\mathbf{x} = (I - A)^{-1}\mathbf{d}.$$

The matrix $L \equiv (I - A)^{-1}$ is called the **Leontief inverse** or **total requirements matrix**. Its $(i,j)$ element $\ell_{ij}$ gives the total output of industry $i$ — direct and indirect — required per unit of final demand for industry $j$'s product. To see why the inverse captures all indirect rounds, expand it as the Neumann series:

$$L = (I - A)^{-1} = I + A + A^2 + A^3 + \cdots$$

The $I$ term represents the direct output; the $A$ term represents first-round intermediate demand; $A^2$ represents the second-round demand triggered by the first round; and so on. The series converges because the spectral radius of $A$ is less than one in a productive economy.

### Output Multipliers

The **column sum** of the Leontief inverse, $\mu_j = \sum_i \ell_{ij}$, is the **output multiplier** of industry $j$: the total output across all industries required to deliver one unit of final demand from industry $j$. Industries with high output multipliers have dense upstream supply chains and trigger large economy-wide output responses when their final demand increases.

The Leontief framework also allows computation of **employment multipliers** (jobs per unit of final demand), **income multipliers** (wages per unit of final demand), and **import multipliers** (foreign inputs induced by domestic final demand). These are constructed by pre-multiplying the Leontief inverse by the appropriate direct coefficient vector:

$$\text{employment multiplier}_j = \mathbf{e}' L \mathbf{e}_j,$$

where $\mathbf{e}$ is the vector of direct labor coefficients (employment per unit of gross output by industry) and $\mathbf{e}_j$ is the $j$th unit vector.

---

## 4.6 Worked Example: A Three-Sector Economy

To make the Leontief framework concrete, consider a stylized three-sector economy: Agriculture (A), Manufacturing (M), and Services (S). The technical coefficient matrix, calibrated loosely to a middle-income economy, is:

$$A = \begin{pmatrix} 0.10 & 0.20 & 0.05 \\ 0.15 & 0.25 & 0.10 \\ 0.10 & 0.15 & 0.20 \end{pmatrix},$$

where entry $a_{ij}$ is the share of industry $j$'s gross output sourced as intermediate inputs from industry $i$. Manufacturing, for example, uses $0.25$ of its own gross output as intermediate inputs and draws $0.20$ from Agriculture.

**Step 1: Form $I - A$.**

$$I - A = \begin{pmatrix} 0.90 & -0.20 & -0.05 \\ -0.15 & 0.75 & -0.10 \\ -0.10 & -0.15 & 0.80 \end{pmatrix}.$$

**Step 2: Compute the Leontief inverse $L = (I - A)^{-1}$** (values rounded to three decimal places):

$$L \approx \begin{pmatrix} 1.168 & 0.367 & 0.113 \\ 0.279 & 1.538 & 0.218 \\ 0.210 & 0.357 & 1.320 \end{pmatrix}.$$

**Step 3: Interpret.** The second column of $L$ says: to deliver \$1 of final demand from Manufacturing, the economy must produce \$0.367 of Agriculture, \$1.538 of Manufacturing (direct plus all indirect rounds), and \$0.357 of Services. The output multiplier for Manufacturing is $\mu_M = 0.367 + 1.538 + 0.357 = 2.262$: each dollar of manufacturing final demand generates \$2.26 of economy-wide gross output.

**Step 4: Policy application.** A public infrastructure procurement program increases Manufacturing final demand by $\Delta d_M = 100$. The required output increase in each sector:

$$\Delta \mathbf{x} = L \cdot (0,\; 100,\; 0)' = (36.7,\; 153.8,\; 35.7)'.$$

Agriculture must expand by \$36.7, Manufacturing by \$153.8, and Services by \$35.7 — even though the direct injection was exclusively into Manufacturing. The input-output framework renders visible the supply-chain interdependencies that an aggregate model obscures entirely.

**Step 5: National accounts consistency check.** Value added per sector is $\mathbf{v} = (I - A)\mathbf{x}$. Then $\sum_i v_i = \mathbf{1}'(I - A)\mathbf{x} = \mathbf{1}'\mathbf{d}$, confirming that total value added equals total final demand — the production approach to GDP equals the expenditure approach, as national accounting requires.

---

## 4.7 National Accounts and Their Limits

The national accounts measure what they were designed to measure — marketed production — with increasing accuracy. But several important limitations deserve careful attention, not to dismiss the accounts but to calibrate how much weight to place on GDP as a welfare measure.

**Non-marketed production.** Household production (childcare, cooking, home repair), volunteer work, and the informal economy are excluded from GDP. In high-income countries, satellite accounts suggest unpaid household production may be equivalent to 25–40% of measured GDP. The exclusion creates systematic biases: when childcare migrates from household to commercial provision, GDP rises even if aggregate social welfare is unchanged.

**Depreciation of natural capital.** Standard accounts subtract the depreciation of manufactured capital $D$ in moving from GNP to NNP. They do not subtract the depletion of natural capital — fisheries, aquifers, soil fertility, or atmospheric carbon absorption capacity. The System of Environmental-Economic Accounting (SEEA), developed alongside the SNA, attempts to fill this gap with satellite natural capital accounts [connects to C:Ch.18]. Countries that grow GDP rapidly by liquidating natural resources appear prosperous by conventional measures while actually running down their productive asset base. Genuine saving — NNP adjusted for natural capital depletion — can be negative even when measured GDP is growing.

**Distribution.** GDP per capita says nothing about how income is distributed. Two economies with the same per capita GDP but Gini coefficients of 0.25 and 0.55 represent radically different social realities. The relationship between macroeconomic aggregates and distributional outcomes is central to Part VIII [Ch. 38].

**Subjective wellbeing.** The Easterlin (1974) paradox documents that beyond a threshold, rising average income is not reliably associated with rising average reported happiness. This does not make GDP useless — it remains the best single predictor of access to education, healthcare, and physical security — but it does mean GDP is a tool with a specific purpose, not a measure of flourishing in any comprehensive sense.

These limitations motivate alternative frameworks: the Human Development Index (HDI), the Genuine Progress Indicator (GPI), Bhutan's Gross National Happiness accounting, and the OECD Better Life Index. None has displaced GDP as the primary summary statistic of macroeconomic performance, partly because of the massive existing infrastructure for GDP measurement and partly because quarterly GDP data provide a uniquely timely and internationally comparable signal of short-run economic conditions for the specific policy questions that macroeconomists most often face.

---

## Chapter Summary

- The circular flow implies the **saving-investment identity** $(S - I) + (T - G) = NX$, an accounting constraint that holds regardless of the state of the economy. It means a current-account deficit must be matched by a private saving deficit, a fiscal deficit, or both.

- The **sectoral financial balances identity** decomposes this into household, corporate, and government sectors. If any sector moves into deficit, at least one other must move toward surplus. This framework provided a precise pre-crisis diagnosis of the U.S. imbalances that contributed to the 2008 crisis.

- The **national income hierarchy** runs from GDP through GNP, NNP, National Income, Personal Income, and Personal Disposable Income, each step stripping out income generated in production that does not reach households' budget constraints. The labor share $\omega = W/Y$ connects this accounting to growth theory and inequality.

- The **balance of payments identity** $CA + FA + KA + \Delta R = 0$ records all cross-border transactions. Current-account deficits require capital inflows or reserve drawdowns, and countries that exhaust reserves without attracting capital face forced adjustment through exchange rate depreciation or a sudden stop.

- The **Leontief input-output model** $\mathbf{x} = (I - A)^{-1}\mathbf{d}$ captures all direct and indirect supply-chain effects of a demand shock, where the Leontief inverse expands as the Neumann series $I + A + A^2 + \cdots$. Column sums give output multipliers; pre-multiplication by factor coefficients gives employment and income multipliers.

- GDP has important limits as a welfare measure: it excludes non-marketed production, ignores natural capital depletion, is insensitive to distribution, and imperfectly proxies subjective wellbeing.

---

## Exercises

**4.1.** In a closed economy with no external sector, household income is $Y = 5{,}000$, consumption $C = 3{,}500$, and government purchases $G = 800$ financed by lump-sum taxes $T = 800$. (a) Compute private saving $S$ and investment $I$. (b) Verify the saving-investment identity. (c) Now the government cuts taxes to $T = 600$ while holding $G = 800$. If households consume 60 cents of every extra dollar of disposable income ($\Delta C = 120$), what must happen to investment if income $Y$ is fixed? Interpret using the sectoral balances identity.

**4.2.** Country A has GDP $= 1{,}000$. Domestic residents earn $80$ in wages and dividends from abroad; foreigners earn $50$ from activities in Country A. Depreciation is $100$. Indirect taxes net of subsidies are $60$. Retained corporate earnings are $40$, corporate taxes $30$, payroll taxes $20$, and transfers received by households $90$. Personal income taxes are $150$. (a) Compute GNP, NNP, National Income, Personal Income, and Personal Disposable Income. (b) If the household saving rate out of PDI is $8\%$, compute household saving in levels.

**4.3.** A country's balance of payments (in billions of dollars): goods exports $= 300$; goods imports $= 360$; net services $= +20$; net primary income $= -15$; net secondary income $= +10$; net FDI inflows $= +25$; net portfolio inflows $= +35$; net other investment outflows $= -20$. (a) Compute the trade balance, the current account, and the financial account. (b) What is the implied change in official reserves? (c) If the central bank wants to prevent currency depreciation, should it buy or sell domestic currency in the foreign exchange market?

**4.4.** ★ Using the three-sector economy of Section 4.6: (a) Final demand shifts to $\mathbf{d} = (200, 150, 400)'$. Compute gross output $\mathbf{x}$. (b) A drought raises the agricultural input requirement for manufacturing from $a_{AM} = 0.20$ to $0.28$, but improves agriculture's own efficiency so $a_{AA}$ falls from $0.10$ to $0.05$. Recompute $A$, $(I - A)$, and the Leontief inverse. Has the output multiplier for manufacturing risen or fallen? (c) Compute the value-added vector $\mathbf{v} = (I - A)\mathbf{x}$ for the original economy at $\mathbf{d} = (200, 300, 500)'$ and verify that $\sum_i v_i = \sum_i d_i$.

**4.5.** ★ The labor share of income is $\omega_A = 0.55$, $\omega_M = 0.40$, $\omega_S = 0.60$, and gross output levels (billions) are $\mathbf{x} = (300, 500, 700)'$. Depreciation is $12\%$ of gross output in each sector. (a) Compute total wage income, total gross operating surplus, NNP, and the net profit share of NNP. (b) The corporate surplus $(S_F - I_F)$ has risen persistently in advanced economies since the 2000s. Using the sectoral financial balances identity and holding the current account constant, what are the implications for the household and government sectors? What macroeconomic conditions would be consistent with this accounting?

**4.6.** ★★ Consider a two-sector economy with sectors $K$ (capital-intensive) and $L$ (labor-intensive) with:

$$A = \begin{pmatrix} 0.30 & 0.10 \\ 0.20 & 0.25 \end{pmatrix}, \qquad \ell_K = 0.25, \quad \ell_L = 0.50$$

where $\ell_j$ is direct labor hours per unit of gross output in sector $j$. (a) Compute the Leontief inverse $L = (I - A)^{-1}$. (b) Compute the full labor content per unit of final demand for each sector, $f_j = \boldsymbol{\ell}' L \mathbf{e}_j$. Which sector has higher full labor content once all indirect supply-chain rounds are counted? (c) Compare the ranking by direct labor intensity ($\ell_K$ vs. $\ell_L$) with the ranking by full labor intensity ($f_K$ vs. $f_L$). Under what structural conditions — captured in $A$ — could a sector with low direct labor intensity have higher full labor intensity than an ostensibly labor-intensive sector? Does this suggest a resolution of the Leontief paradox?

---

*Next: Chapter 5 — Economic Growth and Development: The Long-Run Perspective*
