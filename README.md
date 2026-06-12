# QT

Материалы по соединению **астемизол** (антигистамин, блокатор hERG, удлиняет интервал QT).

## Файлы
- `astemizole.smi` — SMILES структуры.
- `astemizole.md` — идентификаторы (SMILES, InChI, InChIKey, формула, CAS, CID).
- `rxn_export_sequences.py` — выгрузка всех маршрутов ретросинтеза из IBM RXN в CSV.

## Экспорт ретросинтеза (IBM RXN → CSV)

В веб-интерфейсе IBM RXN маршруты («Sequence N») просматриваются по одному, а
единой выгрузки всех маршрутов нет. Скрипт достаёт их все через API в один CSV.

```bash
pip install rxn4chemistry

# Ключ задаётся одним из способов (не коммитится в репозиторий):
export RXN_API_KEY="ваш-ключ"        # либо положите его в локальный файл .rxn_key

python3 rxn_export_sequences.py "COc1ccc(CCN2CCC(Nc3nc4ccccc4n3Cc3ccc(F)cc3)CC2)cc1" --out routes.csv
```

> Файлы `.rxn_key`, `*.env` и `rxn_sequences.csv` исключены через `.gitignore` —
> ключ API никогда не попадёт в git.
