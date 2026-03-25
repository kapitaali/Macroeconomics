# Chapter 4 — The Circular Flow and National Income Accounting

> *"The national income of a country can be looked at in three ways: as the total of incomes earned, as the total of expenditures, or as the total of outputs."*
> — Richard Stone, Nobel Lecture, 1984

---

Chapter 3 introduced the major macroeconomic aggregates in isolation. This chapter asks a more structural question: how do these aggregates fit together? How does money flow through an economy — from firms to households as wages and profits, from households back to firms as consumption, from households and firms to the government as taxes, from the government and abroad back as spending? Understanding these flows is not merely bookkeeping. The accounting identities that govern the circular flow constrain what is macroeconomically possible and reveal connections between the fiscal deficit, private saving, and the current account that are not obvious from any single sector's perspective.

---

## 4.1 The Circular Flow of Income and Expenditure

The simplest macroeconomic model has two sectors: households and firms. Firms hire factor inputs from households — labor and capital — and pay for them as wages and profits. Households use this income to purchase the goods and services firms produce. Money flows in a circle: payments from firms to households as factor income, and payments from households to firms as consumer expenditure.

This circular structure immediately implies a fundamental accounting identity: **aggregate income equals aggregate expenditure**. Every dollar of expenditure is simultaneously someone's income — the payment for a good is the revenue that allows the firm to pay its workers and shareholders. This is not a theory; it is a consequence of double-entry bookkeeping applied at the economy-wide level.

In a more realistic economy with a government and an external sector, the circular flow has **injections** and **leakages**. Injections are flows that add to the demand for domestically produced output beyond what comes from household consumption: investment $I$ (firms spending on new capital), government spending $G$, and exports $X$ (foreign demand for domestic goods). Leakages are flows that remove purchasing power from the domestic circular flow: saving $S$ (income not spent on consumption), taxes $T$, and imports $M$ (spending on foreign rather than domestic goods).

The fundamental circular flow identity for an open economy with a government sector:

$$Y = C + S + T + M = C + I + G + X.$$

Since both expressions equal $Y$, canceling $C$ and rearranging:

$$(S - I) + (T - G) = X - M \equiv NX.$$

This is the **saving-investment identity**. It says that the difference between private saving and private investment, plus the government's fiscal surplus, equals net exports. It holds as an accounting identity regardless of the state of the economy and regardless of any behavioral assumption. Its policy relevance is powerful: a country running persistent current-account deficits ($NX < 0$) must, by accounting necessity, be running either a private saving deficit ($S < I$) or a fiscal deficit ($T < G$), or both. There is no escape from this arithmetic.

---

## 4.2 Sectoral Financial Balances

The saving-investment identity can be extended to a three-sector decomposition that is particularly useful for macroeconomic diagnosis. Denote household sector saving $S_H$, corporate retained earnings $S_F$, government saving (fiscal surplus) $S_G = T - G$, and the current account balance $CA$ (which equals $NX$ plus net factor income from abroad and net transfer payments). The **sectoral financial balances identity** is:

$$\underbrace{(S_H - I_H)}_{\text{household surplus}} + \underbrace{(S_F - I_F)}_{\text{corporate surplus}} + \underbrace{S_G}_{\text{government surplus}} = CA.$$

This identity, developed systematically by Wynne Godley (1999), is a powerful organizing framework for understanding macroeconomic imbalances. Its key implication: if any one sector moves toward deficit, one or more of the remaining sectors must move toward surplus. The government cannot run a larger deficit without either reducing the private sector surplus or worsening the current account, and vice versa. Applied to the U.S. in the 2000s: the simultaneous deterioration of the government surplus and the household surplus (negative $S_H - I_H$, as households borrowed against rising house prices) was arithmetically associated with a large current-account deficit — a global imbalance that contributed to the financial crisis of 2008.

---

## 4.3 The National Income and Product Accounts

The national accounts are the systematic implementation of these identities by statistical agencies. The United States System of National Accounts (SNA 2008) constructs a hierarchy of aggregate income measures, each stripping out a different component:

$$\text{GDP} = \text{GNP} - \text{Net factor income from abroad}$$

$$\text{NNP} = \text{GNP} - \text{Depreciation (D)}$$

$$\text{NI} = \text{NNP} - T^{ind} + \text{Subsidies}$$

$$\text{PI} = \text{NI} - \text{Retained earnings} - \text{Corporate taxes} - \text{Payroll taxes} + \text{Transfers}$$

$$\text{PDI} = \text{PI} - T^{personal}$$

**Definition (Gross National Product).** **GNP** is the total income earned by a country's residents, regardless of where they are located. It differs from GDP — which measures production within a country's borders — by the amount of net factor income from abroad: income that domestic residents earn in other countries minus income that foreign residents earn domestically.

The distinction matters for countries with large numbers of migrant workers (net factor income is positive, so GNP > GDP) or for countries with heavy foreign ownership of domestic capital (net factor income is negative, GNP < GDP). For the United States, GDP and GNP are nearly identical because inflows and outflows of factor income approximately cancel.

**Definition (Net National Product).** **NNP** is GNP minus depreciation — the consumption of fixed capital during the period. NNP is a better measure of the sustainable income a country can produce, because it accounts for the fact that some current production is devoted to replacing worn-out capital rather than adding to the capital stock. GDP is used more commonly in practice because depreciation is harder to measure precisely than gross investment.

**Definition (Personal Disposable Income).** **PDI** is the income available to households after all taxes have been paid and all transfer payments received. It is the relevant budget constraint for households' consumption and saving decisions. The aggregate saving rate $s = S_H / \text{PDI}$ is a central variable in growth theory (Chapter 5) and in debates about the effectiveness of fiscal policy (Chapter 22).

---

## 4.4 The Balance of Payments

The **balance of payments** is the system of accounts that records all economic transactions between residents of one country and residents of all other countries during a given period. It consists of three main accounts, which together sum to zero.

**Definition (Current Account).** The **current account (CA)** records transactions in goods and services (the trade balance), factor income (wages and investment income earned across borders), and transfer payments (remittances, foreign aid). The current account balance is approximately equal to net exports $NX$ for countries with small cross-border factor income flows.

**Definition (Financial Account).** The **financial account (FA)** records cross-border transactions in financial assets: foreign direct investment, portfolio investment (stocks and bonds), and other capital flows. A positive FA means foreigners are accumulating claims on domestic assets (a capital inflow); a negative FA means domestic residents are accumulating claims on foreign assets (a capital outflow).

**Definition (Official Reserve Account).** The **official reserve account ($\Delta R$)** records changes in official foreign exchange reserves held by the central bank. An increase in reserves ($\Delta R > 0$) represents the central bank accumulating foreign assets, typically to defend an exchange rate peg.

The balance of payments identity:

$$CA_t + FA_t + \Delta R_t = 0.$$

This identity holds because every transaction generates two entries — an asset flow and a corresponding payment flow — of equal magnitude. Its policy implication: a current-account deficit must be financed either by net capital inflows ($FA > 0$) or by drawing down reserves ($\Delta R < 0$). Countries that cannot attract sufficient capital inflows while running current-account deficits face balance-of-payments crises — a topic developed in Chapters 26 and 32.

---

## 4.5 Input–Output Analysis

The production-side national accounts aggregate value added across industries. A more disaggregated view is provided by the Leontief (1941) input–output framework, which maps the full web of inter-industry transactions in the economy.

The key concept is the **technical coefficient matrix** $A = [a_{ij}]$, where $a_{ij}$ is the dollar value of input from industry $i$ required to produce one dollar of gross output in industry $j$. Define the gross output vector $\mathbf{x} = (x_1,\ldots,x_n)'$ and the final demand vector $\mathbf{d} = (d_1,\ldots,d_n)'$. Final demand includes private consumption, investment, government spending, and exports. The accounting identity for each industry: gross output equals deliveries to other industries plus deliveries to final demand.

In matrix form: $\mathbf{x} = A\mathbf{x} + \mathbf{d}$. The unique solution (assuming $I - A$ is invertible with non-negative inverse, which holds when the economy is productive):

$$\mathbf{x} = (I - A)^{-1}\mathbf{d}.$$

The matrix $(I-A)^{-1}$ is called the **Leontief inverse**. Its $(i,j)$ element gives the total output of industry $i$ — direct and indirect — required per unit of final demand for industry $j$'s product. "Indirect" here means all the upstream effects: a unit of car demand requires direct steel input, which requires direct iron ore and coal inputs, which require direct energy inputs, and so on.

Input–output multipliers derived from the Leontief inverse allow analysts to trace the economy-wide consequences of sectoral shocks — a supply-chain disruption, a defense spending increase, or a pandemic-induced collapse of tourism — through the full inter-industry network. They provide a structural decomposition of GDP that the single-aggregate models of subsequent chapters cannot offer.

---

*Next: Chapter 5 — Economic Growth and Development: The Long-Run Perspective*
