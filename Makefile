thorns/GaugeWave:	GaugeWave.m
			../Kranc/Bin/kranc GaugeWave.m

clean:
	rm -rf thorns

commit-thorns:	thorns/GaugeWave
		rm -rf master_thorns
		git clone -b master_thorns . master_thorns
		cp -av thorns/* master_thorns/
		git --git-dir $(PWD)/master_thorns/.git --work-tree $(PWD)/master_thorns commit -a -m "Updated from source"
		git fetch master_thorns
