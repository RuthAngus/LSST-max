import numpy as np


def generate_visits(Nvisits=900, tspan=10, stat=False,
                    seasonscale=365./5):
    '''
    Use some very crude approximations for how visits will be spaced out:

    - Survey starts at midnight, time = 0.0
    - Can only observe at night, time > 0.75 | time < 0.25
    - Exposures are clustered around a season w/ a gaussian shape each year
    - Field is observable for first half of year, 0 < date < 182
    - On average, each field should be hit every 3 days during observable season

    Set "stat=True" if you want a plot and a couple basic statistics about the cadence

    '''

    # generate random times for visit, between [0.75 and 0.25]
    time_of_day = np.random.random(Nvisits)/2. - 0.25

    date_of_year = np.floor(np.random.normal(loc=365./4., scale=seasonscale, size=Nvisits))

    year_of_obs = np.floor(np.random.random(Nvisits) * tspan) * 365.

    date_obs = time_of_day + date_of_year + year_of_obs

    date_obs.sort()

    if stat is True:
        print('mean time between visits:')
        print(np.mean(date_obs[1:] - date_obs[:-1]))

        print('median time between visits:')
        print(np.median(date_obs[1:] - date_obs[:-1]))

        plt.figure()
        _ = plt.hist(date_obs, bins=np.arange(date_obs.min(), date_obs.max(),7),
                     histtype='stepfilled', color='k')
        plt.xlabel('Time (days)')
        plt.ylabel('# Visits per Week')
        plt.show()

    return date_obs
